from anytree import Node, RenderTree
import numpy as np
import graphviz as gv
from anytree.dotexport import RenderTreeGraph


class MyNode(Node):
    def __init__(self, name, parent=None):
       Node.__init__(self, name, parent)
       self.name = name
       self.parent=parent
       self.tuple=[]
       self.is_dead=False
    def getTuple(self):
        return self.tuple
    def setDead(self):
        self.is_dead=True

def is_in(a, mytuple):
    if float(a)>float(mytuple[0]) and float(a)<=float(mytuple[1]):
        return True
    return False


n = raw_input("n:")
Beta = raw_input("beta:")
Alpha = raw_input("alpha:")
Delta = raw_input("delta:")



n= int(n) #Number of leaves
Alpha = float(Alpha)
Beta = float(Beta)
Delta = float(Delta)
Ui = np.random.uniform(0.0,1.0,n-1) #Array of n-1 numbers
Vi = np.random.uniform(0.0,1.0,n-1)
Di = np.random.uniform(0.0,1.0,n-1)
Bi = np.random.beta(float(Alpha+1),float(Beta+1),n-1)



Tree = []
Tree.append(MyNode("0: [0,1]"))
Tree[0].tuple=[0,1]
Tree.append(MyNode(str(1)+": [0,"+"{0:.2f}".format(Bi[0])+"]"))
Tree.append(MyNode(str(2)+": ["+"{0:.2f}".format(Bi[0])+",1]"))
Tree[2].parent=Tree[0]
Tree[1].parent=Tree[0]
Tree[1].tuple=[0,Bi[0]]
Tree[2].tuple=[Bi[0],1]

tree_number=2
j=1

while j<n-1:
    if Vi[j] < Delta :
            for tr in Tree:
                if tr.is_leaf and is_in(Di[j], tr.getTuple()):
                    if (not tr.is_dead):
                         tr.name = tr.name+"*"
                    tr.setDead()
                    break
    else:
            for tree in Tree:
                if tree.is_leaf and is_in(Ui[j], tree.getTuple()) and not tree.is_dead:
                    a,b = tree.getTuple()
                    tree_number+=2
                    #Two new children are born here
                    middle = float(Bi[j])*float((float(b)-float(a)))+float(a)
                    Tree.append(MyNode(str(tree_number-1)+": ["+"{0:.4f}".format(a)+","+"{0:.4f}".format(middle)+"]", parent=tree))
                    Tree.append(MyNode(str(tree_number)+":["+"{0:.4f}".format(middle)+","+"{0:.4f}".format(b)+"]", parent=tree))

                    #The new intervals are assigned here
                    Tree[tree_number-1].tuple=[a,middle]
                    Tree[tree_number].tuple=[middle,b]
                    break

    j+=1

for pre, fill, node in RenderTree(Tree[0]):
    print("%s%s" % (pre, node.name))

RenderTreeGraph(Tree[0]).to_picture('beta_splitting.jpg')
