#!/usr/bin/env python3
import numpy as np
import pandas as pd
from pandas.core.frame import DataFrame
from typing import List
from collections import OrderedDict
import random
import sys
import os

import torch.nn as nn
import torch.nn.functional as F
import torch
from torch.utils.data import Dataset
from torch.utils.data import DataLoader
import torch.autograd as autograd

from torch.autograd import Variable
from annoy import AnnoyIndex

disease = sys.argv[1]
#input one disease name (from multiple batches), eg. ACVD, advanced adenoma, CRC, Crohns disease, Healthy, IGT, Obesity, Overweight, Rheumatoid Arthritis, Symptomatic atherosclerosis, T2D, Ulcerative colitis, Underweight

class ScDataset(Dataset):
  def __init__(self):
    self.dataset = []
    self.variable = None
    self.labels = None
    self.transform = None
    self.sample = None
    self.trees = []

  def __len__(self):
    return 10 * 1024

  def __getitem__(self, index):
    dataset_samples = []
    for j, dataset in enumerate(self.dataset):
      rindex1 = np.random.randint(len(dataset))
      rindex2 = np.random.randint(len(dataset))
      alpha = np.random.uniform(0, 1)
      sample = alpha*dataset[rindex1] + (1-alpha)*dataset[rindex2]
      dataset_samples.append(sample)
    return(dataset_samples)

class Encoder(nn.Module):
  def __init__(self, latent_dim):
    super(Encoder, self).__init__()
    self.encoder = nn.Sequential(
      nn.Linear(data_size, 1024),
      nn.BatchNorm1d(1024),
      Mish(),
      nn.Linear(1024, 512),
      nn.BatchNorm1d(512),
      Mish(),
      nn.Linear(512, latent_dim),
    )

  def forward(self, x):
    return(self.encoder(x))

class Decoder(nn.Module):
  def __init__(self, latent_dim):
    super(Decoder, self).__init__()
    self.relu = torch.nn.ReLU()
    self.decoder = nn.Sequential(
      nn.Linear(latent_dim, 512),
      Mish(),
      nn.Linear(512, 1024),
      Mish(),
      nn.Linear(1024, data_size),
    )
    self.decoder2 = nn.Sequential(
      nn.Linear(n_classes, 512),
      Mish(),
      nn.Linear(512, 1024),
      Mish(),
      nn.Linear(1024, data_size),
    )

  def forward(self, ec, es):
    return(self.relu(self.decoder(torch.cat((ec, es), dim=-1))+self.decoder2(es)))

class Mish(nn.Module):
  def __init__(self):
    super().__init__()

  def forward(self, x):
    return (x * torch.tanh(F.softplus(x)))

def weights_init_normal(m):
  classname = m.__class__.__name__
  if(classname.find("Linear") != -1):
    torch.nn.init.normal_(m.weight.data, 0.0, 0.02)
  elif(classname.find("BatchNorm") != -1):
    torch.nn.init.normal_(m.weight.data, 1.0, 0.02)
    torch.nn.init.constant_(m.bias.data, 0.0)

def cat_data(data_A: np.float32, data_B: np.float32, labels: List[List[int]]=None):
  data = np.r_[data_A, data_B]
  if labels is None:
    label = np.zeros(len(data_A)+len(data_B))
    label[-len(data_B):] = 1
    label = np.array([label]).T
  else:
    label = np.r_[labels[0], labels[1]]
  return(data, label)

def setup_seed(seed):
  torch.manual_seed(seed)
  torch.cuda.manual_seed_all(seed)
  np.random.seed(seed)
  torch.backends.cudnn.deterministic = True

key = 'batch'
n_epochs = 10 #150
num_workers=0
lr = 0.0005
b1 = 0.5
b2 = 0.999
latent_dim = 256
n_critic = 5
lambda_co = 3
lambda_rc = 1
seed = 8
setup_seed(seed)

path = 'tables'
train = pd.read_csv('%s/4347_final_relative_abundances.txt' % path, sep='\t', index_col=0).T
metadata = pd.read_csv('%s/Final_metadata_4347.csv' % path, index_col=0).T
train.insert(0,'batch',metadata['Author (year)'].tolist())
train.insert(0,'Type','')
train['Type']=pd.Series(train.index.values).str.split('_',expand=True)[0].tolist()
train_p = train[train['Type']==disease]
batches = sorted(list(set(train['batch'])))
batches_list = [ 'Other' if batch in set(train_p['batch'].value_counts()[6:].index) else batch for batch in train_p['batch'] ]
batches = sorted(list(set(batches_list)))
b = train_p['batch']
train_p['batch'] = batches_list
scd = ScDataset()
scd.variable = np.array(train_p.columns.tolist()[1:])
adata_values = [[train_p[train_p['batch']==batch].iloc[:,2:].values][0] for batch in batches]
std_ = [np.sum(np.std(item, axis=0)) for item in adata_values]
orders = np.argsort(std_)[::-1]
obs_names = [np.array( train_p[train_p['batch']== batch].index.tolist()) for batch in batches]

ec_obs_names = None
for item in orders:
  if ec_obs_names is None:
    ec_obs_names = obs_names[item]
  else:
    ec_obs_names = np.r_[ec_obs_names, obs_names[item]]

scd.dataset = [adata_values[i] for i in orders]
dataloader = DataLoader(
    dataset = scd,
    batch_size=512,
    num_workers=num_workers,
  )

global data_size
global n_classes
data_size = scd.dataset[0].shape[1]
n_classes = len(scd.dataset)

EC = Encoder(latent_dim)
Dec = Decoder(latent_dim + n_classes)
mse_loss = torch.nn.MSELoss()

cuda = True if torch.cuda.is_available() else False
FloatTensor = torch.cuda.FloatTensor if cuda else torch.FloatTensor
LongTensor = torch.cuda.LongTensor if cuda else torch.LongTensor
if cuda:
  EC.cuda()
  Dec.cuda()
  mse_loss.cuda()

EC.apply(weights_init_normal)
Dec.apply(weights_init_normal)
optimizer_Dec = torch.optim.Adam(Dec.parameters(), lr=lr, betas=(b1, b2))
optimizer_EC = torch.optim.Adam(EC.parameters(), lr=lr, betas=(b1, b2))
#Dec.train()

for epoch in range(n_epochs):
  Dec.train()
  EC.train()

  for i, data in enumerate(dataloader):
    datum = [Variable(item.type(FloatTensor)) for item in data]
    batch_size = datum[0].shape[0]

    ES_data1 = -np.zeros((n_classes * batch_size, n_classes))
    for j in range(n_classes):
      ES_data1[j*batch_size:(j+1)*batch_size, j] = 1
    ES_data1 = Variable(torch.tensor(ES_data1).type(FloatTensor))
    ES_data2 = -np.zeros((n_classes * batch_size, n_classes))
    ES_data2[np.arange(n_classes*batch_size),np.random.randint(n_classes, size=n_classes*batch_size)] = 1
    ES_data2 = Variable(torch.tensor(ES_data2).type(FloatTensor))

    optimizer_Dec.zero_grad()
    optimizer_EC.zero_grad()

    loss1_data1 = torch.cat(datum, dim=0)
    loss4 = mse_loss(EC(loss1_data1), EC(Dec(EC(loss1_data1), ES_data2)))
    ae_loss = mse_loss(Dec(EC(loss1_data1), ES_data1), loss1_data1)

    all_loss = (lambda_co * loss4) + (lambda_rc * ae_loss)
    all_loss.backward()

    optimizer_Dec.step()
    optimizer_EC.step()

  print("[Epoch %d/%d] [Reconstruction loss: %f] [Cotent loss: %f]" % (epoch+1, n_epochs, ae_loss.item(), loss4.item(),))

Dec.eval()
EC.eval()

with torch.no_grad():
  data = Variable(FloatTensor(scd.dataset[0]))
  label = np.full((len(scd.dataset[0]),1), batches[orders[0]])
  static_sample = EC(data)
  transform_data = static_sample.cpu().detach().numpy()
  for j in range(1, len(scd.dataset)):
    data = Variable(FloatTensor(scd.dataset[j]))
    static_sample = EC(data)
    fake_data = static_sample.cpu().detach().numpy()
    fake_label = np.full((len(scd.dataset[j]),1), batches[orders[j]])
    transform_data, label = cat_data(transform_data, fake_data, [label, fake_label])

def setup_seed(seed):
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    torch.backends.cudnn.deterministic = True
    random.seed(seed)

def calculate_gradient_penalty(real_data, fake_data, D):
    eta = torch.FloatTensor(real_data.size(0),1).uniform_(0,1)
    eta = eta.expand(real_data.size(0), real_data.size(1))
    cuda = True if torch.cuda.is_available() else False
    if cuda:
        eta = eta.cuda()
    else:
        eta = eta

    interpolated = eta * real_data + ((1 - eta) * fake_data)

    if cuda:
        interpolated = interpolated.cuda()
    else:
        interpolated = interpolated

    # define it to calculate gradient
    interpolated = Variable(interpolated, requires_grad=True)

    # calculate probability of interpolated examples
    prob_interpolated = D(interpolated)

    # calculate gradients of probabilities with respect to examples
    gradients = autograd.grad(outputs=prob_interpolated, inputs=interpolated,
                           grad_outputs=torch.ones(
                               prob_interpolated.size()).cuda() if cuda else torch.ones(
                               prob_interpolated.size()),
                           create_graph=True, retain_graph=True)[0]

    grad_penalty = ((gradients.norm(2, dim=1) - 1) ** 2).mean()
    return grad_penalty

def normalize(data: np.float32) -> np.float32:
    norm = data#(np.exp2(data)-1)
    return norm

def acquire_pairs(X, Y, k, metric):
    X = normalize(X)
    Y = normalize(Y)

    f = X.shape[1]
    t1 = AnnoyIndex(f, metric)
    t2 = AnnoyIndex(f, metric)
    for i in range(len(X)):
        t1.add_item(i, X[i])
    for i in range(len(Y)):
        t2.add_item(i, Y[i])
    t1.build(10)
    t2.build(10)

    mnn_mat = np.bool8(np.zeros((len(X), len(Y))))
    sorted_mat = np.array([t2.get_nns_by_vector(item, k) for item in X])
    for i in range(len(sorted_mat)):
        mnn_mat[i,sorted_mat[i]] = True
    _ = np.bool8(np.zeros((len(X), len(Y))))
    sorted_mat = np.array([t1.get_nns_by_vector(item, k) for item in Y])
    for i in range(len(sorted_mat)):
        _[sorted_mat[i],i] = True
    mnn_mat = np.logical_and(_, mnn_mat)
    pairs = [(x, y) for x, y in zip(*np.where(mnn_mat>0))]
    return pairs

def create_pairs_dict(pairs):
    pairs_dict = {}
    for x,y in pairs:
        if x not in pairs_dict.keys():
            pairs_dict[x] = [y]
        else:
            pairs_dict[x].append(y)
    return pairs_dict

class ScDataset(Dataset):
    def __init__(self, n_sample=3000):
        self.dataset = []
        self.cali_dataset = []
        self.variable = None
        self.anchor_index = 0
        self.query_index = 1
        self.pairs = None
        self.labels = None
        self.transform = None
        self.sample = None
        self.metric = 'euclidean'
        self.k1 = None
        self.k2 = None
        self.n_sample = n_sample


    def change_dataset(self, index: int=1):
        self.query_index = index


    def acquire_anchor(self, index: int=0):
        self.anchor_index = index


    def calculate_mnn_pairs(self):
        tmp = np.arange(len(self.dataset[self.anchor_index]))
        np.random.shuffle(tmp)
        self.sample = self.cali_dataset[self.anchor_index][tmp[:self.n_sample]]
        ####
        tmp2 = np.arange(len(self.dataset[self.query_index]))
        np.random.shuffle(tmp2)
        self.query_sample = self.cali_dataset[self.query_index][tmp2[:self.n_sample]]
        ####
        
        if (self.k1 is None) or (self.k2 is None):
            self.k2 = int(min(len(self.sample), len(self.query_sample))/100)
            self.k1 = max(int(self.k2/2), 1)
        
        print('Calculating Anchor Pairs...')
        anchor_pairs = acquire_pairs(self.sample, self.sample, self.k1, self.metric)
        print('Calculating Query Pairs...')
        query_pairs = acquire_pairs(self.query_sample, self.query_sample, self.k1, self.metric)
        print('Calculating KNN Pairs...')
        pairs = acquire_pairs(self.sample, self.query_sample, self.k1, self.metric)
        print('Calculating Random Walk Pairs...')
        anchor_pairs_dict = create_pairs_dict(anchor_pairs)
        query_pairs_dict = create_pairs_dict(query_pairs)
        pair_plus = []
        for x, y in pairs:
            start = (x, y)
            for i in range(50):
                pair_plus.append(start)
                start = (random.choice(anchor_pairs_dict[start[0]]), random.choice(query_pairs_dict[start[1]]))

        self.datasetA = self.dataset[self.query_index][tmp2[:self.n_sample]][[y for x,y in pair_plus], :]
        self.datasetB = self.dataset[self.anchor_index][tmp[:self.n_sample]][[x for x,y in pair_plus], :]
        print('Done.')

    def __len__(self):
        return 10*1024


    def __getitem__(self, index):
        return random.choice(self.datasetA), random.choice(self.datasetB)


def cat_data(data_A: np.float32, data_B: np.float32, labels: List[List[int]]=None):
    data = np.r_[data_A, data_B]
    if labels is None:
        label = np.zeros(len(data_A)+len(data_B))
        label[-len(data_B):] = 1
        label = np.array([label]).T
    else:
        label = np.r_[labels[0], labels[1]]
    return data, label

def traingan(scd, n_dataset, n_epochs):
    scd.change_dataset(n_dataset)
    scd.calculate_mnn_pairs()

    n_epochs = n_epochs
    n_classes = 2
    data_size = scd.dataset[0].shape[1]
    lr = 0.0002
    b1 = 0.5
    b2 = 0.999
    latent_dim = 256
    n_critic = 100


    cuda = True if torch.cuda.is_available() else False
    FloatTensor = torch.cuda.FloatTensor if cuda else torch.FloatTensor
    LongTensor = torch.cuda.LongTensor if cuda else torch.LongTensor

    dataloader = DataLoader(
        dataset = scd,
        batch_size=1024,
    )

    class Generator(nn.Module):
        def __init__(self):
            super(Generator, self).__init__()
            self.relu = nn.ReLU(inplace=True)
            self.encoder = nn.Sequential(
                nn.Linear(data_size, 1024),
                nn.BatchNorm1d(1024),
                Mish(),
                nn.Linear(1024, 512),
                nn.BatchNorm1d(512),
                Mish(),
                nn.Linear(512, latent_dim),
                nn.BatchNorm1d(latent_dim),
            )
            self.decoder = nn.Sequential(
                nn.Linear(latent_dim, 512),
                nn.BatchNorm1d(512),
                Mish(),
                nn.Linear(512, 1024),
                nn.BatchNorm1d(1024),
                Mish(),
                nn.Linear(1024, data_size),
            )

        def forward(self, x):
            latent_data = self.encoder(x)
            gen_data = self.decoder(latent_data)
            return self.relu(gen_data + x)

    class Discriminator(nn.Module):
        def __init__(self):
            super(Discriminator, self).__init__()
            self.model = nn.Sequential(
                nn.Linear(data_size, 512),
                Mish(),
                nn.Linear(512, 512),
                Mish(),
            )

            # Output layers
            self.adv_layer = nn.Sequential(nn.Linear(512, 1))

        def forward(self, data):
            out = self.model(data)
            validity = self.adv_layer(out)
            return validity

    # Initialize generator and discriminator
    G_AB = Generator()
    D_B = Discriminator()

    if cuda:
        G_AB.cuda()
        D_B.cuda()

    # Initialize weights
    G_AB.apply(weights_init_normal)
    D_B.apply(weights_init_normal)

    optimizer_G_AB = torch.optim.Adam(G_AB.parameters(), lr=lr, betas=(b1, b2))
    optimizer_D_B = torch.optim.Adam(D_B.parameters(), lr=lr, betas=(b1, b2))

    for epoch in range(n_epochs):
        G_AB.train()
        for i, (data_A, data_B) in enumerate(dataloader):
            batch_size = data_A.shape[0]

            # Configure input
            real_data = Variable(data_B.type(FloatTensor))

            # ---------------------
            #  Train Discriminator
            # ---------------------
            optimizer_D_B.zero_grad()
            z = Variable(data_A.type(FloatTensor))
            gen_data = G_AB(z)


            # Loss for real images
            real_validity  = D_B(real_data)
            fake_validity  = D_B(gen_data)


            # Compute W-div gradient penalty
            div_gp = calculate_gradient_penalty(real_data, gen_data, D_B)

            # Adversarial loss
            db_loss = -torch.mean(real_validity) + torch.mean(fake_validity) + 10*div_gp
            db_loss.backward()
            optimizer_D_B.step()

            # -----------------
            #  Train Generator
            # -----------------

            if i % n_critic == 0:
                optimizer_G_AB.zero_grad()
                z = Variable(data_A.type(FloatTensor), requires_grad=True)
                gen_data = G_AB(z)
                fake_validity = D_B(gen_data)
                gab_loss = -torch.mean(fake_validity)
                gab_loss.backward()

                optimizer_G_AB.step()


        # --------------
        # Log Progress
        # --------------

        print(
            "[Epoch %d/%d] [D loss: %f] [G loss: %f]"
            % (epoch+1, n_epochs,
               db_loss.item(),
               gab_loss.item(),
              )
        )

    G_AB.eval()
    with torch.no_grad():
        z = Variable(FloatTensor(scd.dataset[scd.query_index]))
        static_sample = G_AB(z)
        fake_data = static_sample.cpu().detach().numpy()
    return fake_data


cali_batches = sorted(list(set(label.T[0])))
scd1 = scd
scd = ScDataset(len(label.T[0]))
scd.metric = 'angular'
scd.k1 = None
scd.k2 = None
scd.variable = scd1.variable

print('Orders:','<-'.join(batches[i] for i in orders))
scd.dataset = [adata_values[i] for i in orders]
cali_adata = pd.DataFrame(transform_data)
cali_adata.insert(0,'batch',label.T[0])
cali_adata_values = [np.array(cali_adata[cali_adata['batch'] == batch].iloc[:,1:]) for batch in batches]
cali_orders = orders
scd.cali_dataset = [cali_adata_values[i] for i in cali_orders]
scd.transform = np.copy(scd.dataset[scd.anchor_index])

for i in range(1, len(scd.dataset)):
  print(f'Merging dataset {batches[orders[i]]} to {batches[orders[0]]}')
  fake_data = traingan(scd, i, n_epochs=n_epochs)
  scd.transform = np.r_[scd.transform, fake_data]

output = pd.DataFrame(scd.transform)
output.columns=train_p.columns.tolist()[2:]
out_index =[]
for i in orders:
  out_index += train_p[train_p['batch']==batches[i]].index.tolist() 
output.insert(0,'ID',out_index)

metadata.insert(0,'ID',train.index.tolist())
output.insert(1,'batch',output['ID'].replace(train['batch'].to_dict()))
output.insert(1,'Type',output['ID'].str.split('_',expand=True)[0])
output.to_csv('demo/output/%s-BER.csv' % disease) #save output file, eg. Overweight-BER.csv
print('done!')
os._exit(0)

