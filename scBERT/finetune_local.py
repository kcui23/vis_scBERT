# -*- coding: utf-8 -*-
import os
import argparse
import random
import torch
from torch import nn
from torch.optim import Adam
from torch.utils.data import DataLoader, Dataset
import scanpy as sc
import numpy as np
import pickle as pkl
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.metrics import accuracy_score, f1_score, confusion_matrix, classification_report

from performer_pytorch import PerformerLM
from utils import *  # Make sure this includes necessary functions like seed_all

os.chdir("./scBERT-1.0.0")

parser = argparse.ArgumentParser()
parser.add_argument("--bin_num", type=int, default=4, help='Number of bins.')
parser.add_argument("--gene_num", type=int, default=16906, help='Number of genes.')
parser.add_argument("--epoch", type=int, default=100, help='Number of epochs.')
parser.add_argument("--seed", type=int, default=2021, help='Random seed.')
parser.add_argument("--batch_size", type=int, default=3, help='Batch size.')
parser.add_argument("--learning_rate", type=float, default=1e-4, help='Learning rate.')
parser.add_argument("--valid_every", type=int, default=1, help='Validation frequency.')
parser.add_argument("--pos_embed", type=bool, default=True, help='Use Gene2vec encoding.')
parser.add_argument("--data_path", type=str, default='./data/Zheng68k.h5ad', help='Data path.')
parser.add_argument("--model_path", type=str, default='./model/panglao_pretrained.pth', help='Pretrained model path.')
parser.add_argument("--ckpt_dir", type=str, default='./ckpts/', help='Checkpoint directory.')
parser.add_argument("--model_name", type=str, default='finetuned', help='Finetuned model name.')

args = parser.parse_args()

SEED = args.seed
EPOCHS = args.epoch
BATCH_SIZE = args.batch_size
LEARNING_RATE = args.learning_rate
SEQ_LEN = args.gene_num + 1

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
seed_all(SEED)

class SCDataset(Dataset):
    def __init__(self, data, label):
        super().__init__()
        self.data = data
        self.label = label

    def __getitem__(self, index):
        rand_start = random.randint(0, self.data.shape[0]-1)
        full_seq = self.data[rand_start].toarray()[0]
        full_seq[full_seq > (args.bin_num + 1)] = args.bin_num + 1
        full_seq = torch.from_numpy(full_seq).long()
        full_seq = torch.cat((full_seq, torch.tensor([0]))).to(device)
        seq_label = self.label[rand_start]
        return full_seq, seq_label

    def __len__(self):
        return self.data.shape[0]

class Identity(torch.nn.Module):
    def __init__(self, dropout = 0., h_dim = 100, out_dim = 10):
        super(Identity, self).__init__()
        self.conv1 = nn.Conv2d(1, 1, (1, 200))
        self.act = nn.ReLU()
        self.fc1 = nn.Linear(in_features=SEQ_LEN, out_features=512, bias=True)
        self.act1 = nn.ReLU()
        self.dropout1 = nn.Dropout(dropout)
        self.fc2 = nn.Linear(in_features=512, out_features=h_dim, bias=True)
        self.act2 = nn.ReLU()
        self.dropout2 = nn.Dropout(dropout)
        self.fc3 = nn.Linear(in_features=h_dim, out_features=out_dim, bias=True)

    def forward(self, x):
        x = x[:,None,:,:]
        x = self.conv1(x)
        x = self.act(x)
        x = x.view(x.shape[0],-1)
        x = self.fc1(x)
        x = self.act1(x)
        x = self.dropout1(x)
        x = self.fc2(x)
        x = self.act2(x)
        x = self.dropout2(x)
        x = self.fc3(x)
        return x

data = sc.read_h5ad(args.data_path)
label_dict, label = np.unique(np.array(data.obs['celltype']), return_inverse=True)
with open('label_dict.pkl', 'wb') as fp:
    pkl.dump(label_dict, fp)
with open('label.pkl', 'wb') as fp:
    pkl.dump(label, fp)
label = torch.from_numpy(label)
data = data.X

sss = StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=SEED)
for train_index, test_index in sss.split(np.zeros(data.shape[0]), label.numpy()):
    data_train, data_val = data[train_index], data[test_index]
    label_train, label_val = label[train_index], label[test_index]

train_dataset = SCDataset(data_train, label_train)
val_dataset = SCDataset(data_val, label_val)

train_loader = DataLoader(train_dataset, batch_size=BATCH_SIZE, shuffle=True)
val_loader = DataLoader(val_dataset, batch_size=BATCH_SIZE)

model = PerformerLM(
    num_tokens=args.bin_num + 3,
    dim=200,
    depth=6,
    max_seq_len=SEQ_LEN,
    heads=10,
    local_attn_heads=0,
    g2v_position_emb=args.pos_embed
)

model_path = args.model_path
ckpt = torch.load(model_path, map_location=device)
model.load_state_dict(ckpt['model_state_dict'])
model.to_out = Identity(dropout=0., h_dim=128, out_dim=label_dict.shape[0])
model = model.to(device)

optimizer = Adam(model.parameters(), lr=LEARNING_RATE)

for epoch in range(EPOCHS):
    model.train()
    for batch in train_loader:
        data, labels = batch
        optimizer.zero_grad()
        outputs = model(data)
        loss = nn.CrossEntropyLoss()(outputs, labels)
        loss.backward()
        optimizer.step()

    if epoch % args.valid_every == 0:
        model.eval()
        with torch.no_grad():
            correct = 0
            total = 0
            for data, labels in val_loader:
                outputs = model(data)
                _, predicted = torch.max(outputs.data, 1)
                total += labels.size(0)
                correct += (predicted == labels).sum().item()
        print('Accuracy of the network on validation data: %d %%' % (100 * correct / total))

# Save your model
torch.save(model.state_dict(), os.path.join(args.ckpt_dir, args.model_name + '.pth'))

