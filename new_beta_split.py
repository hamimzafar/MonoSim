from trees import *
import numpy as np
import random

def is_in(a, mytuple):
    if float(a)>float(mytuple[0]) and float(a)<=float(mytuple[1]):
        return True
    return False


n = raw_input("n:")
Beta = raw_input("beta:")
Alpha = raw_input("alpha:")
Delta = raw_input("delta:")

n = int(n)  # Number of leaves
Alpha = float(Alpha)
Beta = float(Beta)
Delta = float(Delta)
Ui = np.random.uniform(0.0,1.0,n-1)  # Array of n-1 numbers
Vi = np.random.uniform(0.0,1.0,n-1)
Di = np.random.uniform(0.0,1.0,n-1)
Bi = np.random.beta(float(Alpha+1), float(Beta+1), n-1)

Tree = []
Tree.append(FullBiTree("0: [0,1]"))
#Tree[0].tuple=[0,1]
Tree[0].set_node_property('tuple', [0,1])
Tree.append(FullBiTree(str(1)+": [0,"+"{0:.2f}".format(Bi[0])+"]"))
Tree.append(FullBiTree(str(2)+": ["+"{0:.2f}".format(Bi[0])+",1]"))
Tree[1].set_node_property('tuple',[0,Bi[0]])
Tree[2].set_node_property('tuple',[Bi[0],1])
Tree[0].set_node_property('dead',False)
Tree[1].set_node_property('tuple',False)
Tree[2].set_node_property('tuple',False)
Tree[0].set_children(Tree[1], Tree[2])
for i in range(3):
    rand = random.random
    Tree[i].set_node_property('branch length',rand)

tree_number=2
j=1

while j<n-1:
    if Vi[j] < Delta :
            for tr in Tree:
                if tr.is_leaf and is_in(Di[j], tr.get_node_property('tuple')):
                    if not tr.get_node_property('dead'):
                        tr.name = tr.name+"*"
                    tr.set_node_property('dead',True)
                    break
    else:
            for tree in Tree:
                if tree.is_leaf and is_in(Ui[j], tree.get_node_property('tuple')) and not tree.get_node_property('dead'):
                    a,b = tree.get_node_property('tuple')
                    tree_number += 2
                    #Two new children are born here
                    middle = float(Bi[j])*float((float(b)-float(a)))+float(a)
                    lchild = FullBiTree(str(tree_number-1)+": ["+"{0:.4f}".format(a)+","+"{0:.4f}".format(middle)+"]")
                    rchild = FullBiTree(str(tree_number) + ":[" + "{0:.4f}".format(middle) + "," + "{0:.4f}".format(b) + "]")
                    lchild.set_node_property('branch length', random.random)
                    rchild.set_node_property('branch length', random.random)
                    tree.set_children(lchild, rchild)
                    Tree.append(FullBiTree(str(tree_number-1)+": ["+"{0:.4f}".format(a)+","+"{0:.4f}".format(middle)+"]"))
                    Tree.append(FullBiTree(str(tree_number)+":["+"{0:.4f}".format(middle)+","+"{0:.4f}".format(b)+"]"))

                    #The new intervals are assigned here
                    Tree[tree_number-1].set_node_property('tuple', [a,middle])
                    Tree[tree_number].set_node_property('tuple',[middle,b])
                    break

    j+=1

