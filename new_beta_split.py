from trees import *
import numpy as np
import random

def is_in(a, mytuple):
    if float(a)>float(mytuple[0]) and float(a)<=float(mytuple[1]):
        return True
    return False


# n = raw_input("n:")
# Beta = raw_input("beta:")
# Alpha = raw_input("alpha:")
# Delta = raw_input("delta:")

def beta_split(n, Beta, Alpha, Delta):
    Ui = np.random.uniform(0.0,1.0,n-1)  # Array of n-1 numbers
    Vi = np.random.uniform(0.0,1.0,n-1)
    Di = np.random.uniform(0.0,1.0,n-1)
    Bi = np.random.beta(float(Alpha+1), float(Beta+1), n-1)

    Treelist = []
    Treelist.append(FullBiTree("0: [0,1]"))
    #Tree[0].tuple=[0,1]
    Treelist[0].set_node_property('tuple', [0,1])
    Treelist.append(FullBiTree(str(1)+": [0,"+"{0:.2f}".format(Bi[0])+"]"))
    Treelist.append(FullBiTree(str(2)+": ["+"{0:.2f}".format(Bi[0])+",1]"))
    Treelist[1].set_node_property('tuple',[0,Bi[0]])
    Treelist[2].set_node_property('tuple',[Bi[0],1])
    Treelist[0].set_node_property('leaf',False)
    Treelist[1].set_node_property('leaf',True)
    Treelist[2].set_node_property('leaf',True)
    #Treelist[0].set_node_property('dead', False)
    #Treelist[1].set_node_property('dead', False)
    #Treelist[2].set_node_property('dead', False)
    Treelist[0].set_children(Treelist[1], Treelist[2])
    for i in range(3):
        rand = random.random
        Treelist[i].set_node_property('branch length',rand)

    tree_number=2
    j=1


    while j<n-1:
        #if Vi[j] < Delta:
        #        for tr in Treelist:
        #            if tr.is_leaf and is_in(Di[j], tr.get_node_property('tuple')):
        #                if not tr.get_node_property('dead'):
        #                    tr.name = tr.name+"*"
        #                tr.set_node_property('dead',True)
        #                break
        #else:
        for tree in Treelist:
            #print Tree[1]
            if tree.get_node_property('leaf') and is_in(Ui[j], tree.get_node_property('tuple')):
                #and not tree.get_node_property('dead'):
                a,b = tree.get_node_property('tuple')
                tree_number += 2
                #Two new children are born here
                middle = float(Bi[j])*float((float(b)-float(a)))+float(a)
                lchild = FullBiTree(str(tree_number-1)+": ["+"{0:.4f}".format(a)+","+"{0:.4f}".format(middle)+"]")
                rchild = FullBiTree(str(tree_number) + ":[" + "{0:.4f}".format(middle) + "," + "{0:.4f}".format(b) + "]")
                lchild.set_node_property('branch length', random.random)
                rchild.set_node_property('branch length', random.random)
                lchild.set_node_property('leaf', True)
                rchild.set_node_property('leaf', True)
                tree.set_node_property('leaf', False)
                tree.set_children(lchild, rchild)
                #print tree
                Treelist.append(lchild)
                Treelist.append(rchild)

                #The new intervals are assigned here
                Treelist[tree_number - 1].set_node_property('tuple', [a,middle])
                Treelist[tree_number].set_node_property('tuple',[middle,b])
                break
        j += 1
    return Treelist[0]

