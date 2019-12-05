'''
Created on 2017年10月14日

@author: Song
'''
import py2neo
import math
import string
#from py2neo.types import  Subgraph
#from py2neo import Node,Graph,Relationship, NodeSelector
from py2neo import Graph

class myneo4jconnect():
    def connectGraph(self):
        py2neo.authenticate("localhost:7474", USERNAME, PASSWORD)
        graph = Graph("http://localhost:7474/db/data")
        return graph   
  
    def CreateRelationship(self):     #Interaction   
        graph = self.connectGraph()      
        #for i in range(1,2):
            #s="%s-ADGenome.txt" %(i) 
        for i in range(1,6):  
            file = open(r'..\chr%d_4DG-hg38.txt' % (i))
            for line in file:
                line=line[:-1]          #去掉换行符
                #print(line)
                a=line.split('\t')           #根据'\t'分割字符串 
                cypher_Creat ="""MATCH (n1:%s_Range {Start:'%s',End:'%s'}) MATCH (n2:%s_Range {Start:'%s',End:'%s'}) create  (n1)-[:Interaction {CellType:'%s',Method:'%s',Confidence:'%s',PubMedID:'%s',SourceDB:'CistromeDB'}]->(n2)"""  % ( a[0],a[1].strip(),a[2].strip(),a[3],a[4].strip(),a[5].strip(),a[6].strip(),a[7].strip(),a[8].strip(),a[9].strip())        
                
                #cypher_Creat ="""MATCH (n1:%s_Range {Start:'%s',End:'%s'}) -[r:Interaction {CellType:'%s'}]- (n2:%s_Range {Start:'%s',End:'%s'}) set r.Confidence='%s'"""  % ( a[0],a[1].strip(),a[2].strip(),a[6].strip(),a[3],a[4].strip(),a[5].strip(),a[8].strip())        
                
                graph.run(cypher_Creat)
                #print(cypher_Creat)   
                #break     
            file.close()
            #break
             
    def CreateRangeIndex(self):        
        graph = self.connectGraph()
        for num in range(1,23):
            cypher_Creat ="""CREATE INDEX ON :chr%d_Range(Start)"""% num              
            graph.run(cypher_Creat)
        cypher_Creat ="""CREATE INDEX ON :chrX_Range(Start)"""          
        graph.run(cypher_Creat)
        cypher_Creat ="""CREATE INDEX ON :chrY_Range(Start)"""            
        graph.run(cypher_Creat) 
    
    def CreateBinIndex(self):        
        graph = self.connectGraph()
        for num in range(1,23):
            cypher_Creat ="""CREATE INDEX ON :chr%d(Name)"""% num              
            graph.run(cypher_Creat)
        cypher_Creat ="""CREATE INDEX ON :chrX(Name)"""          
        graph.run(cypher_Creat)
        cypher_Creat ="""CREATE INDEX ON :chrY(Name)"""            
        graph.run(cypher_Creat) 
   
    def CreateBRRelationship(self): #Includsion
        graph = self.connectGraph()      
        #for i in range(1,2):
        #a=1    #s="%s-ADGenome.txt" %(i)   
        for num in range(1,23):   #1，2，3，4，5\13-22
            file = open(r'..\%d-gene-hg38.csv'% (num))
            for line in file:            
                line=line[:-1]          #去掉换行符
                
                a=line.split(',')      
                for i in range(math.floor(int(a[0])/200000+1),math.floor(int(a[1])/200000+1)+1):
                    cypher_Creat ="""MATCH (n1:chr%d {BinName:'Bin%d'}) MATCH (n2:chr%d_Range {Start:'%s',End:'%s'}) create  (n1)-[:Inclusion]->(n2)"""  % (num,i,num,a[0].strip(),a[1].strip())        
                    graph.run(cypher_Creat)
                    #print(cypher_Creat)
            #break
        
    def DeleteIndex(self):        
        graph = self.connectGraph()
        for num in range(1,23):
            cypher_Creat ="""drop index on :chr%d(Start)"""% num              
            graph.run(cypher_Creat)
        cypher_Creat ="""drop index on :chrX(Start)"""          
        graph.run(cypher_Creat)
        cypher_Creat ="""drop index on :chrY(Start)"""            
        graph.run(cypher_Creat) 
        
    def CreateRangeNodes(self):        
        graph = self.connectGraph()
        for num in range(1,23):
            cypher_Creat ="""LOAD CSV WITH HEADERS FROM "file:///%d-gene-hg38.csv" AS row CREATE (n:chr%d_Range) SET n = row"""% (num,num)              
            graph.run(cypher_Creat)
        cypher_Creat ="""LOAD CSV WITH HEADERS FROM "file:///X-gene-hg38.csv" AS row CREATE (n:chrX_Range) SET n = row"""
        graph.run(cypher_Creat) 
        cypher_Creat ="""LOAD CSV WITH HEADERS FROM "file:///Y-gene-hg38.csv" AS row CREATE (n:chrY_Range) SET n = row"""     
        graph.run(cypher_Creat)           
    
    def CreateBinNodes(self):        
        graph = self.connectGraph()
        for num in range(1,23):
            cypher_Creat ="""LOAD CSV WITH HEADERS FROM "file:///h%d-ADGenome-hg38.csv" AS row CREATE (n:chr%d) SET n = row"""% (num,num)              
            graph.run(cypher_Creat)
        cypher_Creat ="""LOAD CSV WITH HEADERS FROM "file:///hX-ADGenome-hg38.csv" AS row CREATE (n:chrX) SET n = row"""
        graph.run(cypher_Creat) 
        cypher_Creat ="""LOAD CSV WITH HEADERS FROM "file:///hY-ADGenome-hg38.csv" AS row CREATE (n:chrY) SET n = row"""     
        graph.run(cypher_Creat)
    
    def PrintRelationship(self):        
        file = open(r'..\1-ADGenome.txt')
        for line in file:            
            line=line[:-1]          #去掉换行符
            a=line.split('\t')           #根据'\t'分割字符串
            cypher_Creat ="""MATCH (n1:%s {BinName:'%s'}) MATCH (n2:%s {BinName:'%s'}) create (n1)-[:Interaction {Start:'%s',End:'%s',CellType:'%s',PubMedID:'%s',SourceDB:'CistromeDB'}]->(n2)"""             
            graph = self.connectGraph()
            graph.run(cypher_Creat)
        
    def AddProtienRelationship(self):        
        graph = self.connectGraph()       
        #for i in range(1,2):
            #s="%s-ADGenome.txt" %(i)        
        file = open(r'..\h3.txt')
        for line in file:            
            line=line[:-1]          #去掉换行符
            a=line.split('\t')           #根据'\t'分割字符串
            file2 = open(r'..\%s' % (a[7]))
            
            #cypher_Creat ="""create (n:Protein {CellType:'%s',Name:'%s',GEO:'%s'})"""% (a[3],a[4],a[5],a[6],a[1])              
            #graph.run(cypher_Creat)
            for line2 in file2:
                line2=line2[:-1]          #去掉换行符
                b=line2.split('\t')
                if len(b[0])<=5:
                    cypher_Creat ="""MATCH (n1:%s {Bins:'Bin%d'}) MATCH (n2:Protein {Name:'%s',GEO:'%s'}) create (n2)-[:Bind{Start:'%s',End:'%s',SourceDB:'CistromeDB'}]->(n1)"""% (b[0],math.ceil(int(b[1])/2000),a[6],a[1],b[1],b[2])              
                    graph.run(cypher_Creat)
            file2.close()                        
        file.close()
                   
    def AddProtienRelationship2(self):        #Bind
        graph = self.connectGraph()       
        
        for x in range(1,25):
            file2 = open(r'..\%dK562.txt' % x )              
            
            for line2 in file2:
                line2=line2[:-1]          #去掉换行符
                b=line2.split('\t')
                #print(c)
                #c=c+1
                t = math.ceil(int(b[1])/2000)
                x = math.ceil(int(b[2])/2000)
                    #MATCH (n1:%s{CellType=~'.*HELA*,Bins:'Bin%d'}) where n1.CellType=~'.*HELA*.'
               
                for i in range(t,x+1):                                              
                    #cypher_Creat ="""MATCH (n1:%s{Name:'Bin%d'})  MATCH (n2:TF {Name:'%s'}) create (n2)-[:Bind{CellType:'A549',Start:'%s',End:'%s',GEO:'%s',OtherGEO:'%s',Confidence:'Qvalue:%s',SourceDB:'CistromeDB'}]->(n1)"""   % (b[0],i,b[5],b[1].strip(),b[2].strip(),b[4],b[7],b[6].strip())          
                    
                    cypher_Creat ="""MATCH (n1:%s{Name:'Bin%d'})  MATCH (n2:TF {Name:'%s'}) create (n2)-[:Bind{CellType:'%s',Start:'%s',End:'%s',GEO:'%s',OtherGEO:'%s',Confidence:'Qvalue:%s',SourceDB:'CistromeDB'}]->(n1)"""   % (b[0],i,b[5],b[3].strip(),b[1].strip(),b[2].strip(),b[4],b[7],b[6].strip())          
                   
                    #cypher_Creat =""" MATCH (n2:TF {Name:'%s'})-[r:Bind{CellType:'%s'}]->(n1:%s{Name:'Bin%d'}) set r.Confidence='qvalue:%s'"""   % (b[5],b[3].strip(),b[0],i,b[6].strip())          
                    #cypher_Creat ="""CALL apoc.periodic.iterate('MATCH (n:LinkedinID) WHERE n.userDefinedImageUrl IS NOT NULL RETURN n','WITH {n} AS n SET n:头像',{batchSize:10,parallel:true})"""
                    
                    #print(cypher_Creat)
                    graph.run(cypher_Creat)  
                    
                #break                                  
            file2.close() 
            #break
        
         
    def AddTFRelationship(self):  #TFInteraction
        graph = self.connectGraph()      
        
        file2 = open(r'..\TF-TF-comb.csv' )                 
        for line2 in file2:
            #line2=line2[:-1]          #去掉换行符
            b=line2.split(',')
            #print(b)
                                      
            cypher_Creat ="""MATCH (n1:TF{Name:'%s'})  MATCH (n2:TF {Name:'%s'}) create (n1)-[:Interaction{Information:'%s'}]->(n2)"""   % (b[0].strip(),b[1].strip(),b[2])          
            #cypher_Creat ="""match(n:TF{Name:'%s'})-[r:Interaction]-(m:TF{Name:'%s'}) set r.Confidence='%s' """  % (b[2],b[3],b[6]) 
            graph.run(cypher_Creat)  
            #print(cypher_Creat)
            #break                  
        file2.close()
    
    def AddBinRelationship(self):  
        graph = self.connectGraph()      
        a=[90738,85374,79618,72530,69145,66891,67538,66621,57175,53441,50991,45113,41622,40125,29301,32165,23347,25402]                
        for j in range(1,2):
            for i in range(1,28609): #  28609
                cypher_Creat ="""MATCH (n1:chrY{BinName:'Bin%d'})  MATCH(n2:chrY{BinName:'Bin%d'}) merge (n1)-[:Connection]-(n2)"""   % (i,i+1)          
                graph.run(cypher_Creat)  
                #print(cypher_Creat)
                #break                  
            #break      
    def UpdateBRRelationship(self): 
        graph = self.connectGraph()      
        #for i in range(1,2):
            #s="%s-ADGenome.txt" %(i)   

        file = open(r'..\All-ADGenome-hg38.txt')
        i=1
        for line in file:            
            line=line[:-1] 
            i=i+1         #去掉换行符
                #print(line)
            a=line.split('\t')      
            cypher_Creat ="""MATCH (n1:%s_Range {Start:'%s',End:'%s'}) -[r:Interaction{CellType:'%s'}]-(n2:%s_Range {Start:'%s',End:'%s'}) set r.method='%s'"""  % (a[0],a[1].strip(),a[2].strip(),a[6].strip(),a[3],a[4].strip(),a[5].strip(),a[7])        
            graph.run(cypher_Creat) 
            print(i)    
            #break 
    
    def DeleteTFRelationship(self):  
        graph = self.connectGraph()      
        
        file2 = open(r'..\delete-tf.txt' )                 
        for line2 in file2:
            line2=line2[:-1]          #去掉换行符
            b=line2.split('\t')
            #print(b[0])
                                      
            cypher_Creat ="""MATCH (n:TF{Name:'%s'})-[r:TFInteraction]- (m:TF) 
            CREATE (m)-[r2:TFInteraction]->(m) 
            SET r2 = r
            WITH r,n
            DELETE r,n """   % (b[0])          
            graph.run(cypher_Creat)  
            #print(cypher_Creat)
            #break                  
        file2.close()
        
    def DeleteBindRelationship(self):  
        graph = self.connectGraph()      
        
        file2 = open(r'..\TF.csv' )                 
        for line2 in file2:
            line2=line2[:-1]          #去掉换行符
            b=line2.split('\t')
            #print(b[0])
            for i in range(1,23):
                cypher_Creat ="""MATCH (n:TF{Name:'%s'})-[r:Bind]- (m:chrX)  DELETE r """   % (b[0],i)          
                graph.run(cypher_Creat)  
                #print(cypher_Creat)
            #break                  
        file2.close()
    
    def AddRNAProtienRelationship(self):        #Bind
        graph = self.connectGraph()       
        
        file2 = open(r'..\LR_high_hg38-all.txt' )                 
        c=1    
        for line2 in file2:
            line2=line2[:-1]          #去掉换行符
            b=line2.split('\t')
            print(c)
            c=c+1
            t = math.ceil(int(b[8])/2000)
            x = math.ceil(int(b[9])/2000)
                    #MATCH (n1:%s{CellType=~'.*HELA*,Bins:'Bin%d'}) where n1.CellType=~'.*HELA*.'
               
            for i in range(t,x+1): 
                pid = b[7] + ':' + b[8].strip() + '-' + b[9].strip()                                             
                cypher_Creat ="""MATCH (n1:%s{Name:'Bin%d'})  MATCH (n2:LncRNA {Name:'%s'}) create (n2)-[:Bind{LncRNA_ID:'%s',CellType:'%s',Associated_Factors:'%s',Epigenetic_Modifications:'%s',Genomic_region:'%s',PubMedID:'%s',Method:'%s',Confidence:'High-throughput',SourceDB:'LnChrom'}]->(n1)"""   % (b[7],i,b[1].strip(),b[0].strip(),b[5].strip(),b[2].strip(),b[3].strip(),pid,b[6].strip(),b[4].strip())          
                graph.run(cypher_Creat)  
                #print(cypher_Creat)
                #break                                  
        file2.close() 
      
              
a = myneo4jconnect()


#a.CreateBinNodes()
#a.CreateRangeNodes()
#a.CreateBinIndex()
#a.CreateRangeIndex()
#a.CreateRelationship()     #Interaction 
#a.CreateBRRelationship()   #Inclusion

#a.AddProtienRelationship2()  #TF bind
#a.AddTFRelationship()       # TF Interaction
#a.AddRNAProtienRelationship() #LncRNA bind
