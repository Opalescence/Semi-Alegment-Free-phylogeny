#!/usr/bin/env pythonw
import glob
import math
import os
import numpy as np
import pickle
import subprocess
import screed
import shutil
import operator as op
import time
from joblib import Parallel, delayed
import multiprocessing
import matplotlib.pyplot as plt
import numpy as np
from itertools import islice
from sklearn import mixture
import io
import tempfile
import random
import progressbar
from time import sleep

########## LOCATION ###########

#dir="29Escherichia_Exact/1x"
#location='/Users/Phanucheep/data/NGS2/'+dir
########## Functions ##########

def nCr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)

def tofasta(dir):
	datafile = glob.glob('/Users/Phanucheep/data/NGS/'+dir+'/*.fna')
	
	if(len(glob.glob('/Users/Phanucheep/data/NGS/'+dir+'/*.fasta'))!=0):
		return
	print("Make fasta.....")
	for file in datafile:
		i=0;
		name=os.path.splitext(os.path.basename(file))[0]
		si,non=name.split('-')
		sa='/Users/Phanucheep/data/NGS/'+dir+'/'+si+'.fasta'
		query = open(sa,"w+")
	
		with screed.open(file) as seqfile:
			#query.write(">"+name+" "+str(i+1)+".\n")
			for record in seqfile:
				s=record.sequence
				if i == 1000000:
					break
				query.write(">"+name+" "+str(i+1)+".\n")
				query.write(s+"\n\n")
				#query.write(s)
				i=i+1
		
		print(sa+" "+si) #c2compute
		#print(si)	#CV list
	#for file in datafile:
	#	name=os.path.splitext(os.path.basename(file))[0]
	#	si,non=name.split('-')
	#	print(si)
	
	print("Finish")

def to_Database(dataset):
	if(len(glob.glob(location+'/*.udb'))!=0):
		return
	for file in dataset:
		name=os.path.splitext(os.path.basename(file))[0]
		si,non=name.split('-')
		print("Building database for "+si+"...........")
		subprocess.call(["/Users/Phanucheep/./usearch","-makeudb_usearch",file,"-output",location+"/"+si+".udb"])
def count(dataset):
	total=0
	total_len=0
	for file in dataset:
		with screed.open(file) as seqfile:
			for record in seqfile:
				total_len+=len(record.sequence)
				total+=1
	print("Read num :"+str(total))
	print("Total Len :"+str(total_len))
	print()

	
def to_query(dataset):
	if(len(glob.glob(location+'/original_query.fasta'))!=0):
		return
	n=0
	N=len(dataset)
	query=[]
	for file in dataset:
		print(file)
		if(file!="query" and file!="original_query"):
			with screed.open(file) as seqfile:
				i=0
				for record in seqfile:
					query.insert(n+(i*(n+1)),record)	
					i=i+1
			n=n+1
	sa=location+'/original_query.fasta'
	query_file = open(sa,"w+")
	i=0
	for record in query:
		query_file.write(">"+record.name+"\n")
		query_file.write(record.sequence+"\n\n")
		i=i+1
		#if(i>len(query)/2):
		#	break
		
def file_parse(search_results):
	#search_results = glob.glob(location+"/result_"+type+"_*.m8")
	result_num = [[0 for i in range(len(search_results))] for i in range(len(search_results))]
	all_result=[]
	for result in search_results:
		with open(result,'r') as f:
			result_list=[]
			for x in f:
				result_item=x.strip().split('\t') 
				result_list.append(result_item)
			all_result.append(result_list)
			
	return all_result

def get_specie_list(dataset):
	specie_list={}
	specie=0
	for db in dataset:
		tmp=os.path.splitext(os.path.basename(db))[0]
		name=tmp.strip().split('-') 
		if(name[0]!="query" and name[0]!="original_query"):
			specie_list[name[0]]=specie
			specie+=1

	print(specie_list)
	return specie_list
	

def outdist(D,search_results):
	print(str(len(D))+"x"+str(len(D[0])))

	file = open(location+"/distV6_final","w+")

	print("Writing....len")
	file.write("	"+str(len(D))+"\n")

	for i in list(range(len(search_results))):
		tmp=os.path.splitext(os.path.basename(search_results[i]))[0]
		name=tmp.strip().split('_') 
		st=name[2]+"    "
		#print(str(i)+" "+s)
		for di in D[i]:
			st=st+"  "+str(di)
		file.write(st+"\n")
		print(st) 	
		
def outdist_2(D,dataset):
	print(str(len(D))+"x"+str(len(D[0])))

	file = open(location+"/final_dist_test_3","w+")

	print("Writing....len")
	file.write("	"+str(len(D))+"\n")

	for i in list(range(len(dataset))):
		tmp=os.path.splitext(os.path.basename(dataset[i]))[0]
		name=tmp.strip().split('.') 
		st=name[0]+"    "
		#print(str(i)+" "+s)
		for di in D[i]:
			st=st+"  "+str(di)
		file.write(st+"\n")
		print(st) 		
			
def dist_cal_V1(dataset):
	specie_list=get_specie_list(dataset)
	
	if(len(glob.glob(location+"/result_"+type+"_*.m8"))==0):
		test_search(query,type)
		
	search_results = glob.glob(location+"/result_"+type+"_*.m8")
	all_result=file_parse(search_results)

	result_num = [[0 for i in range(len(search_results))] for i in range(len(search_results))]
	total_al_col = [[0 for i in range(len(search_results))] for i in range(len(search_results))]

	d= [[0 for i in range(len(search_results))] for i in range(len(search_results))]
	D= [[0 for i in range(len(search_results))] for i in range(len(search_results))]

	for result in all_result:
		for result_item in result:		
			q_names,q_num=result_item[0].split(' ')
			tmp=q_names.strip().split('-')
			q_name=tmp[0]
			t_names,t_num=result_item[1].split(' ')
			tmp=t_names.strip().split('-')
			t_name=tmp[0]
			
			id=result_item[2]
			al_col=result_item[3]
			q_len=result_item[4]
			t_len=result_item[5]
			
			result_num[specie_list[q_name]][specie_list[t_name]]+=1
			total_al_col[specie_list[q_name]][specie_list[t_name]]+=int(al_col)

	print(result_num)
	print(str(len(result_num))+"x"+str(len(result_num[0])))
	print(total_al_col)
	print(str(len(total_al_col))+"x"+str(len(total_al_col[0])))

	for result in all_result:
		for result_item in result:		
			q_names,q_num=result_item[0].split(' ')
			tmp=q_names.strip().split('-')
			q_name=tmp[0]
			t_names,t_num=result_item[1].split(' ')
			tmp=t_names.strip().split('-')
			t_name=tmp[0]
			
			id=float(result_item[2])
			al_col=int(result_item[3])
			#q_len=result_item[4]
			#t_len=result_item[5]
			
			d[specie_list[q_name]][specie_list[t_name]]+=((1-(id/100))*(al_col/total_al_col[specie_list[q_name]][specie_list[t_name]]))

	for i in list(range(len(search_results))):
		for j in list(range(len(search_results))):
			if(i==j):
				D[i][j]=0
			else:
				total=result_num[i][j]+result_num[j][i]
				D[i][j]=(d[i][j]*result_num[i][j]/total)+(d[j][i]*result_num[j][i]/total)
	
	outdist(D,search_results)

def search(query,type):
	
	shutil.copyfile(query,location+'/query.fasta')  
	query=location+'/query.fasta'
	
	for db in dataset:
			
		name=os.path.splitext(os.path.basename(db))[0]
		
		if(name!="query" and name!="original_query"):
			#subprocess.call(["/Users/Phanucheep/./usearch","-usearch_"+type,query,"-strand","plus","-db",db,"-id","0.25","-query_cov","0.5","-maxaccepts","2","-userout",location+"/result_"+type+"_"+name+".m8","-userfields","query+target+id+alnlen+ql+tl+qlo+qhi+tlo+thi"])
			subprocess.call(["/Users/Phanucheep/./usearch","-usearch_"+type,query,"-strand","plus","-db",db,"-id","0.25","-query_cov","0.5","-userout",location+"/result_"+type+"_"+name+".m8","-userfields","query+target+id+alnlen+ql+tl+qlo+qhi+tlo+thi"])
			#update_query3(location+"/result_"+type+"_"+name+".m8",query)

def search_parallel(query,type,db):
	name=os.path.splitext(os.path.basename(db))[0]
	if(name!="query" and name!="original_query"):
		subprocess.call(["/Users/Phanucheep/./usearch","-usearch_"+type,query,"-strand","plus","-db",db,"-id","0.25","-query_cov","0.5","-userout",location+"/result_"+type+"_"+name+".m8","-userfields","query+target+id+alnlen+ql+tl+qlo+qhi+tlo+thi"])
				
def new_search(query,type):
	
	shutil.copyfile(query,location+'/query.fasta')  
	query=location+'/query.fasta'
	
	Parallel(n_jobs=8)(delayed(search_parallel)(query,type,dataset[i]) for i in range(len(dataset)))		


def histrogram(item):

	print(np.std(item))
	print(np.mean(item))
	n, bins, patches = plt.hist(item, bins=100, color='#0504aa',alpha=0.7, rwidth=0.85)
	plt.grid(axis='y', alpha=0.75)
	plt.xlabel('Similarity of alignment pair')
	plt.ylabel('Frequency')
	plt.title('')
	plt.show()
	

	
def search_parallel_2(query,type,db,index,i):
	
	if(i!=0):
		with open(query) as myfile:
			head = [next(myfile) for x in range(index[i-1])]
		s=""
		for line in head:
			s=s+line
		tmp = tempfile.NamedTemporaryFile()
		with open(tmp.name, 'w') as f:
			f.write(s)
	
		name=os.path.splitext(os.path.basename(db[i]))[0]
		if(name!="query" and name!="original_query"):
			subprocess.call(["/Users/Phanucheep/./usearch","-usearch_"+type,tmp.name,"-strand","plus","-db",db[i],"-id","0.25","-query_cov","0.5","-userout",location+"/result_"+type+"_"+name+".m8","-userfields","query+target+id+mid+pairs+gaps+alnlen+mism+ids"])
			#subprocess.call(["/Users/Phanucheep/./usearch","-usearch_"+type,tmp.name,"-strand","plus","-db",db[i],"-id","0.25","-query_cov","0.5","-userout",location+"/result_"+type+"_"+name+".m8","-userfields","query+target+id+alnlen+ql+tl+qlo+qhi+tlo+thi"])

def cvTree_data(dataset):
	index=[]
	i=0
	list=""
	infile=""
	for file in dataset:
		print(file)
		name=os.path.splitext(os.path.basename(file))[0] 
		sa=location+'/'+name+'.ffn'
		list=list+name+"\n"
		infile=infile+name+".ffn\n"
		ffn_file = open(sa,"w+")
		with screed.open(file) as seqfile:
			for record in seqfile:
				ffn_file.write(">"+record.name+"\n")
				ffn_file.write(record.sequence+"\n\n")

	sa=location+'/list'
	list_file=open(sa,"w+")
	list_file.write(list)
		
	sa=location+'/infile'
	in_file=open(sa,"w+")
	in_file.write(infile)	
		
				
def query_search(dataset):
	if(len(glob.glob(location+'/query.fasta'))!=0):
		return
	N=len(dataset)
	query=[]
	index=[]
	i=0
	for file in dataset:
		print(file)
		if(file!="query" and file!="original_query"):
			with screed.open(file) as seqfile:
				for record in seqfile:
					query.append(record)	
					i=i+1
		
		index.append(3*i)
		
	sa=location+'/query.fasta'
	query_file = open(sa,"w+")
	for record in query:
		query_file.write(">"+record.name+"\n")
		query_file.write(record.sequence+"\n\n")
	
	
	
	query=location+'/query.fasta'
	Parallel(n_jobs=8)(delayed(search_parallel_2)(query,type,dataset,index,i) for i in range(len(dataset)))		
	

def dist_cal_V7(dataset):
	specie_list=get_specie_list(dataset)
	
	if(len(glob.glob(location+"/result_"+type+"_*.m8"))==0):
		query_search(dataset)
	
	search_results = glob.glob(location+"/result_"+type+"_*.m8")
	all_result=file_parse(search_results)

	result_num = [[0 for i in range(len(dataset))] for i in range(len(dataset))]
	total_w = [[0 for i in range(len(dataset))] for i in range(len(dataset))]

	d= [[0 for i in range(len(dataset))] for i in range(len(dataset))]
	D= [[0 for i in range(len(dataset))] for i in range(len(dataset))]
	
	item=[]
	
	for result in all_result:
		for result_item in result:		
			id=result_item[2]
			if(float(id)!=100):
				item.append(float(id)/100)
	item=np.array(item)
	
	histrogram(item)

	T=0.0015/np.var(item)
	#T=0.1
	i=0
	for result in all_result:
		for result_item in result:	
			i+=1	
			q_names=result_item[0].strip().split(' ')[0]
			tmp=q_names.strip().split('-')
			tmp=tmp[0].strip().split('.')
			q_name=tmp[0]
			t_names=result_item[1].strip().split(' ')[0]
			tmp=t_names.strip().split('-')
			tmp=tmp[0].strip().split('.')
			t_name=tmp[0]
			
			id=result_item[2]
			
			result_num[specie_list[q_name]][specie_list[t_name]]+=1

			x=((float(id)/100)-1)/T
			total_w[specie_list[q_name]][specie_list[t_name]]+=math.exp(x)
	

	weight_dis=[]
	
	for result in all_result:
		for result_item in result:	
			q_names=result_item[0].strip().split(' ')[0]
			tmp=q_names.strip().split('-')
			tmp=tmp[0].strip().split('.')
			q_name=tmp[0]
			t_names=result_item[1].strip().split(' ')[0]
			tmp=t_names.strip().split('-')
			tmp=tmp[0].strip().split('.')
			t_name=tmp[0]
			
			id=result_item[2]
			
			dis=jukescantor(float(id))

			x=((float(id)/100)-1)/T
			
			w=math.exp(x)/total_w[specie_list[q_name]][specie_list[t_name]]
			weight_dis.append(w)
				
			d2=dis*w
			
			d[specie_list[q_name]][specie_list[t_name]]+=d2
	#histrogram(item)
	#histrogram(weight_dis)
	for i in list(range(len(dataset))):
		for j in list(range(len(dataset))):
			if(i==j):
				D[i][j]=0
			else:
				total=result_num[i][j]+result_num[j][i]
				D[i][j]=(d[i][j]*result_num[i][j]/total)+(d[j][i]*result_num[j][i]/total)
	outdist_2(D,dataset)


###########################################################################################

def dist_cal_V8(dataset):
	specie_list=get_specie_list(dataset)
	
	if(len(glob.glob(location+"/result_"+type+"_*.m8"))==0):
		query_search()
	
	search_results = glob.glob(location+"/result_"+type+"_*.m8")
	all_result=file_parse(search_results)

	result_num = [[0 for i in range(len(dataset))] for i in range(len(dataset))]
	total_al_col = [[0 for i in range(len(dataset))] for i in range(len(dataset))]
	total_sim = [[0 for i in range(len(dataset))] for i in range(len(dataset))]
	total_w = [[0 for i in range(len(dataset))] for i in range(len(dataset))]
	total_w2 = [[0 for i in range(len(dataset))] for i in range(len(dataset))]

	d= [[0 for i in range(len(dataset))] for i in range(len(dataset))]
	D= [[0 for i in range(len(dataset))] for i in range(len(dataset))]
	
	item=[]
	
	for result in all_result:
		for result_item in result:		
			id=result_item[2]
			if(float(id)!=100):
				item.append(float(id)/100)
	item=np.array(item)
	item=item.reshape(-1, 1)
	
	
	clf = mixture.GaussianMixture(n_components=2, covariance_type='full')
	clf.fit(item)
	cov=clf.covariances_.tolist()
	test=clf.predict_proba(np.array(0.9).reshape(-1, 1))
	p=test.tolist()[0][0]
	
	if(p>0.5):
		q=0
	else:
		q=1
	i=0
	cov=math.sqrt(cov[q][0][0])
	print(cov)
	print("Finish Train......")

	#histrogram(item)
	id=np.arange(0,1,0.01)
	y=[]
	for i in id:
		test=clf.predict_proba(np.array(i).reshape(-1, 1))
		if (i<0.4):
			y.append(0)
		else:
			tmp=test.tolist()[0][q]
			y.append(tmp*math.exp(i))
	
	plt.plot(id,y)
	plt.show()
	T=0.0015/np.var(item)
	T=cov*2
	prob_list=[]
	for result in all_result:
		print(i)
		i+=1
		for result_item in result:			
			q_names=result_item[0].strip().split(' ')[0]
			tmp=q_names.strip().split('-')
			tmp=tmp[0].strip().split('.')
			q_name=tmp[0]
			t_names=result_item[1].strip().split(' ')[0]
			tmp=t_names.strip().split('-')
			tmp=tmp[0].strip().split('.')
			t_name=tmp[0]
			
			id=result_item[2]
			al_col=result_item[3]
			q_len=result_item[4]
			t_len=result_item[5]
			
			result_num[specie_list[q_name]][specie_list[t_name]]+=1
			total_al_col[specie_list[q_name]][specie_list[t_name]]+=int(al_col)
			total_sim[specie_list[q_name]][specie_list[t_name]]+=float(id)

			x=((float(id)/100)-1)/T
			prob=clf.predict_proba(np.array((float(id)/100)).reshape(-1, 1))
			prob=prob.tolist()[0][q]
			prob_list.append(prob)
			total_w[specie_list[q_name]][specie_list[t_name]]+=(math.exp(x)*prob)

	

	weight_dis=[]
	d_dis=[]
	i=0
	for result in all_result:
		for result_item in result:		
			q_names=result_item[0].strip().split(' ')[0]
			tmp=q_names.strip().split('-')
			tmp=tmp[0].strip().split('.')
			q_name=tmp[0]
			t_names=result_item[1].strip().split(' ')[0]
			tmp=t_names.strip().split('-')
			tmp=tmp[0].strip().split('.')
			t_name=tmp[0]
			
			id=result_item[2]
			al_col=result_item[3]
			
			dis=jukescantor(float(id))
			#dis=1-(float(id)/100)
			prob=prob_list[i]
			i+=1
			x=((float(id)/100)-1)/T

			w=(math.exp(x)*prob)/total_w[specie_list[q_name]][specie_list[t_name]]
			
		
			weight_dis.append(w)
			d2=dis*w
			d_dis.append(d2)
			d[specie_list[q_name]][specie_list[t_name]]+=d2

	for i in list(range(len(dataset))):
		for j in list(range(len(dataset))):
			if(i==j):
				D[i][j]=0
			else:
				total=result_num[i][j]+result_num[j][i]
				D[i][j]=(d[i][j]*result_num[i][j]/total)+(d[j][i]*result_num[j][i]/total)
	outdist_2(D,dataset)



def jukescantor(id):
	p=1-(id/100)
	if p<0.75:
		dis=(-3/4)*math.log(1-(4/3*p))
	else:
		dis=p
	dis=p
	return dis


def weight_f_1(x):
  return x**math.e


def dist_cal_Vtest(dataset):
	specie_list=get_specie_list(dataset)
	
	if(len(glob.glob(location+"/result2_"+type+"_*.m8"))==0):
		query_search()
	
	search_results = glob.glob(location+"/result2_"+type+"_*.m8")
	
	all_result=file_parse(search_results)

	result_num = [[0 for i in range(len(dataset))] for i in range(len(dataset))]
	total_w = [[0 for i in range(len(dataset))] for i in range(len(dataset))]

	d= [[0 for i in range(len(dataset))] for i in range(len(dataset))]
	D= [[0 for i in range(len(dataset))] for i in range(len(dataset))]
	
	item=[]
	
	for result in all_result:
		for result_item in result:		
			id=result_item[3]
			
			#di=((float(result_item[5]))+(2*float(result_item[7])))/(float(result_item[6])+float(result_item[7]))
			#id=(1-di)*100
			
			if(float(id)!=100):
				item.append(float(id)/100)
	item=np.array(item)

	T=0.0015/np.var(item)
	#T=0.1
	i=0
	for result in all_result:
		for result_item in result:	
			i+=1	
			q_names=result_item[0].strip().split(' ')[0]
			tmp=q_names.strip().split('-')
			tmp=tmp[0].strip().split('.')
			q_name=tmp[0]
			t_names=result_item[1].strip().split(' ')[0]
			tmp=t_names.strip().split('-')
			tmp=tmp[0].strip().split('.')
			t_name=tmp[0]
			
			id=result_item[3]
			
			#di=((float(result_item[5]))+(2*float(result_item[7])))/(float(result_item[6])+float(result_item[7]))
			#id=(1-di)*100
			
			result_num[specie_list[q_name]][specie_list[t_name]]+=1

			x=((float(id)/100)-1)/T
			total_w[specie_list[q_name]][specie_list[t_name]]+=math.exp(x)
	

	weight_dis=[]
	
	for result in all_result:
		for result_item in result:	
			q_names=result_item[0].strip().split(' ')[0]
			tmp=q_names.strip().split('-')
			tmp=tmp[0].strip().split('.')
			q_name=tmp[0]
			t_names=result_item[1].strip().split(' ')[0]
			tmp=t_names.strip().split('-')
			tmp=tmp[0].strip().split('.')
			t_name=tmp[0]
			
			id=result_item[3]
			
			#di=((float(result_item[5]))+(2*float(result_item[7])))/(float(result_item[6])+float(result_item[7]))
			#id=(1-di)*100
			
			dis=jukescantor(float(id))

			x=((float(id)/100)-1)/T
			
			w=math.exp(x)/total_w[specie_list[q_name]][specie_list[t_name]]
			weight_dis.append(w)
				
			d2=dis*w
			
			d[specie_list[q_name]][specie_list[t_name]]+=d2
	#histrogram(item)
	#histrogram(weight_dis)
	for i in list(range(len(dataset))):
		for j in list(range(len(dataset))):
			if(i==j):
				D[i][j]=0
			else:
				total=result_num[i][j]+result_num[j][i]
				D[i][j]=(d[i][j]*result_num[i][j]/total)+(d[j][i]*result_num[j][i]/total)
	outdist_2(D,dataset)




########## Main ##########

type="global"
dset="0.001"
#dset=""
start_time = time.time()

dir=dset+"_illumina"
location='/Users/Phanucheep/data/NGS2/'+dir
dataset = glob.glob(location+'/*.fq')
start_time = time.time()
print(dataset)

print("Make Query.....")
#query_search(dataset)
print("Finish query....")
#cvTree_data(dataset)
dist_cal_V7(dataset)

