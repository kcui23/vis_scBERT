import requests
from bs4 import BeautifulSoup
from tqdm import tqdm
import pickle

def get_gene_link(url):
    response = requests.get(url)
    webpage = response.content

    soup = BeautifulSoup(webpage, 'html.parser')

    for link in soup.find_all('a'):
        if 'ath:A' in str(link):
            return link

def get_NT_seq(url):
    '''
    return: name, NT seq
    '''
    response = requests.get(url)
    webpage = response.content

    soup = BeautifulSoup(webpage, 'html.parser')
    a=soup.find_all('tr')
    start = str(a[5]).find('cel">')
    end = str(a[5]).find('<br')
    name = str(a[5])[start+5:end]
    if ',' in name:
        name = name.split(',')[0]
    tmp=a[-1]
    tmp=str(tmp).split('\n')
    tmp_dna=''
    for i in tmp[2:]:
        tmp_dna+=i.replace('<br/>','').replace('</td></tr>', '')
    return name, tmp_dna.upper()

def pathway_to_NT_seq(pathway_gene_name_filename, gene_link_K_page_filename, gene_nt_seq_link_filename, nt_seq_dic_pkl_filename):
    with open(pathway_gene_name_filename, 'r') as f:
        gene_ids = []
        for i in f.readlines():
            gene_ids.append(i[:6])

    urls = []
    prefix = 'https://www.genome.jp/entry/'
    for i in gene_ids:
        urls.append(prefix+i)

    gene_links = []
    for i in tqdm(urls, desc='Getting gene links'):
        a = get_gene_link(i)
        start = str(a).find('"')
        end = str(a)[start+1:].find('"')
        a_link = 'https://www.genome.jp'+str(a)[start+1:start+end+1]
        gene_links.append(a_link)

    nt_seq = {}
    for i in tqdm(gene_links, desc='Getting NT sequences'):
        name, tmp_dna = get_NT_seq(i)
        if len(str(name)) < 10 and str(name)!='non':
            nt_seq[name] = tmp_dna
    
    with open(gene_link_K_page_filename, 'wb') as pickle_file:
        pickle.dump(urls, pickle_file)
    
    with open(gene_nt_seq_link_filename, 'wb') as pickle_file:
        pickle.dump(gene_links, pickle_file)
    
    with open(nt_seq_dic_pkl_filename, 'wb') as pickle_file:
        pickle.dump(nt_seq, pickle_file)