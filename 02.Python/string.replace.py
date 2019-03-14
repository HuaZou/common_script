#%%
class String_Replace:
    def __init__(self,mystr,mydic):
        self.str = mystr
        self.dic = mydic
        self.res = []
        self.tmp = []
    def Replace_str(self):
        self.tmp = []
        for i in self.res:
            for j in i:
                if j in self.dic:
                    tmp = []
                    for x in self.dic[j]:
                        tmp.append(i.replace(j,x,1))
                    self.tmp+=(tmp)
                    break
        self.res = self.tmp
    def run(self):
        self.res.append(self.str)
        for myid in self.str:
            if myid in self.dic:
                self.Replace_str()
        return(self.res)
                
#%%
a = "AGCTAGCTYATMA"
b = {"Y":["a","b"],"M":["c","d"]}
c = String_Replace(a,b)
c.run()