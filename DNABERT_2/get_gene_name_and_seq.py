# TODO: 
# many fucntions are duplicated, need to clean up, can be generalized and merged into one function

import random
import requests
from bs4 import BeautifulSoup
from tqdm import tqdm
import re
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


def pathway2NTseq_hsa(pathway_list, output_name='hsa_NT_seq_dict', save=True):
    '''
    Input: 
    - pathway_list, list of pathway names
    - output_name, str, output file name
    - save, bool, save the output or not
    
    Output:
    - list of dict, dict: {gene_name: NT_seq}
    '''
    def get_hsa_gene_link(url):
        l=[]
        r = requests.get(url)
        soup = BeautifulSoup(r.content, 'html.parser')
        for i in soup.find_all('a'):
            href = i.get('href')
            if href is not None and href.startswith('/entry/hsa:'):
                l.append(href)
        return l
    
    prefix = "https://www.kegg.jp/entry/pathway+"
    path_dict = {}
    for i in tqdm(range(len(pathway_list)), desc='Getting gene links'):
        link_i = prefix + pathway_list[i]
        path_dict[pathway_list[i]] = get_hsa_gene_link(link_i)
    
    result = {}
    for k,v in tqdm(path_dict.items(), desc='Getting NT sequences'):
        nt_dict = {}
        for i in v:
            nt_link = 'https://www.kegg.jp'+i
            name = i.split('/')[-1]
            _, seq = get_NT_seq(nt_link)
            nt_dict[name] = seq
        result[k] = nt_dict
    
    if save:
        with open(output_name+'.pkl', 'wb') as f:
            pickle.dump(result, f)

def pathway2NTseq_cmiu(pathway_list, output_name='cmiu_NT_seq_dict', save=True, lmt_len_mode=True):
    '''
    Input: 
    - pathway_list, list of pathway names
    - output_name, str, output file name
    - save, bool, save the output or not
    
    Output:
    - list of dict, dict: {gene_name: NT_seq}
    '''
    def get_cmiu_gene_link(url):
        l=[]
        r = requests.get(url)
        soup = BeautifulSoup(r.content, 'html.parser')
        for i in soup.find_all('a'):
            href = i.get('href')
            if href is not None and href.startswith('/entry/cmiu:'):
                l.append(href)
        return l
    
    prefix = "https://www.kegg.jp/entry/pathway+"
    path_dict = {}
    for i in tqdm(range(len(pathway_list)), desc='Getting gene links'):
        link_i = prefix + pathway_list[i]
        l = get_cmiu_gene_link(link_i)
        if lmt_len_mode and len(l) > 50:
            continue
        path_dict[pathway_list[i]] = l
    
    result = {}
    for k,v in tqdm(path_dict.items(), desc='Getting NT sequences'):
        nt_dict = {}
        for i in v:
            nt_link = 'https://www.kegg.jp'+i
            name = i.split('/')[-1]
            _, seq = get_NT_seq(nt_link)
            nt_dict[name] = seq
        result[k] = nt_dict
    
    if save:
        with open(output_name+'.pkl', 'wb') as f:
            pickle.dump(result, f)
    
    return result

def clean_raw_pathway_str(raw, name=''):
    l=[]
    for i in raw.split('\n'):
        if i:
            p = i.split()[0]
            if p.startswith(name):
                p = p[len(name):]
            p = name + p
            l.append(p)
    return l

def get_tidy_brite(brite_name):
    '''
    *Input*: brite_name: str, e.g. cmiu04131 
    *Output*: soup1: str, soup2: dict
    '''
    url = "https://rest.kegg.jp/get/br:" + brite_name
    r = requests.get(url)
    soup1 = BeautifulSoup(r.text, 'html.parser')
    soup1 = soup1.text
    r = requests.get(url+'/json')
    soup2 = r.json()
    return soup1 ,soup2

def get_nt_seq_tidy(gene_name, prefix='cmiu:'):
    url = "https://rest.kegg.jp/get/" + prefix + gene_name
    r = requests.get(url)
    s = BeautifulSoup(r.text, 'html.parser').text
    split_soup = s.split('\n')
    for i in range(len(split_soup)):
        if split_soup[i].startswith('NT'):
            raw_s = ''.join(split_soup[i+1:])
            return(re.sub(r"[ /]", "", raw_s)).upper()
    
    return None

def find_neighbors(gene_name:str, brite_str:str):
    '''
    **Inputs**:
    - gene_name: str, e.g. B1H56_02580
    - brite_str: str, e.g. brite_str, _ = get_tidy_brite('cmiu00001')
    
    **Outputs**:
    - neighbors: dict, keys are the hierarchies, values are the gene names. 
    hier be like '09100 Metabolism -> 09101 Carbohydrate metabolism -> 00500 Starch and sucrose metabolism [PATH:cmiu00500] -> B1H56_01825 carbohydrate kinase\tK00847 E2.7.1.4; fructokinase [EC:2.7.1.4]'
    '''
    neighbors = {}
    level = ['A','B']
    hier = ''
    names = []
    
    for i in brite_str.split('\n'):
        if i:
            tmp = i[1:].strip().split()
            
            flag = False # flag to indicate whether it is time to update the hier
            if i[0] in ('A', 'B', 'C', 'D', 'E', 'F', 'G'):
                level[1] = i[0]
                flag = True
            if level[1] == 'A':
                hier = i[1:].strip()
            elif flag:
                if level[1]>level[0]:
                    hier += ' -> ' + i[1:].strip()
                    if tmp and tmp[0].startswith(gene_name[:3]):
                        names.append(tmp[0])
                if level[1]==level[0]:
                    # back to the previous level
                    hier = ''.join([i+' -> ' for i in hier.split(' -> ')[:-1]]) + i[1:].strip()
                    if tmp and tmp[0].startswith(gene_name[:3]):
                        names.append(tmp[0])
                if level[1]<level[0]:
                    # back to the previous level twice
                    if gene_name in names:
                        neighbors[hier] = names
                    hier = ''.join([i+' -> ' for i in hier.split(' -> ')[:-2]]) + i[1:].strip()
                    names = []
            level[0] = level[1]
    
    return neighbors

def build_sentences(gene_name, neighbors, prefix='cmiu:', sentence_length=None):
    
    def _single_sentence_group(gene_name, local_neighbor, prefix):
        sent = []
        indices = []
        n_sentences = _get_n_sentences(local_neighbor, sentence_length)
        
        for t in range(n_sentences):
            if len(local_neighbor) < 5:
                names = random.sample(local_neighbor, len(local_neighbor))
            else:
                if len(local_neighbor) < 50:
                    vol = int(len(local_neighbor) * .6)
                    while True:
                        names = random.sample(local_neighbor, vol)
                        if gene_name in names:
                            break
            nt_seqs = []
            for i in names:
                nt_seq = get_nt_seq_tidy(i, prefix)
                if nt_seq:
                    nt_seqs.append(nt_seq)
                    indices.append(local_neighbor.index(i))
            tmp_sentence = ' '.join(nt_seqs)
            sent.append(tmp_sentence)
        return sent, indices

    def _get_n_sentences(local_neighbor, sentence_length):
        # TODO: implement a better way to determine the number of sentences
        if len(local_neighbor) < 5:
            return 1
        return 5

    sentences = {}
    for k, v in tqdm(neighbors.items()):
        sent, indices = _single_sentence_group(gene_name, v, prefix)
        sentences[k] = (sent, indices)
    return sentences