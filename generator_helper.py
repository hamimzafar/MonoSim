import numpy as np
from copy import deepcopy
import copy
from collections import defaultdict
import random
from trees import *

def define_nummutes(t, total_length,numsites):
    """
    :param t: A phylogenetic tree with branch lengths already assigned
    :param total_length: The total branch length of the aforementioned tree
    :return: The same tree with all of teh branch length normalized
    """
    if t.is_leaf() == False:
        bl = t.get_node_property('branch length')
        # the branch length of the branch from the parent of this node to this node
        nbl = float(float(bl)/float(total_length))
        # the normalized branch length to be used to find the number of mutations
        nm = round(numsites * nbl)
        # the number of mutations that took place from the parent of this node to this node
        t.set_node_property('nummutes', nm)
        # setting the number of mutations for this node as a node property
        define_nummutes(t.get_left_child(), total_length, numsites)
        define_nummutes(t.get_right_child(), total_length, numsites)
        # recursively call the function on the rest of the tree
    else:
        bl = t.get_node_property('branch length')
        nbl = float(float(bl) / float(total_length))
        nm = round(numsites * nbl)
        t.set_node_property('nummutes', nm)
        # same as the if statement but it does not call the function recursively because we are at a list

def set_alleles(t, changed_sites, parent):
    """
    :param lt: the left child of the tree that we are about to add this child to
    :param rt: the right child of the tree that we are about to add this child to
    :return: the modified children with a new sequence of 0's and 1's representing their alleles
    """
    mutated_sites = []
    # the spots of the alleles that have been mutated, is used as an input for the recursive call so that we know
    # how many sites have been mutated throughout the tree
    flag = False
    # when this flag turns true, all of the sites that can be mutated have already been mutated
    sites = changed_sites
    # the sites that have been changed so far in the tree
    original_one = parent.get_node_property('copy 1 alleles')
    original_two = parent.get_node_property('copy 2 alleles')
    original = parent.get_node_property('alleles')
    # the allele string of the parent of this node
    one = original_one
    two = original_two
    new = original
    nm = t.get_node_property('nummutes')
    # the predetermined numebr of mutations that will occur on this node
    for num in range(int(nm)):
        # this loop mutates "nm" sites in the string
        count = 0
        site = random.randint(0, len(one) - 1)
        #the site to be mutated
        while site in sites:
            #repeat this process until the site being mutated has not already been mutated because of infinite alleles
            #assumption
            count +=1
            site = random.randint(0, len(one) - 1)
            if len(sites) == len(original_one):
                flag = True
                break
        if flag == False:
            sites.append(site)
            mutated_sites.append(site)
            # add the site to the mutated sites of the tree so that we know for the rest of the nodes
            rand2 = random.random()
            if rand2 > .5:
                two = two[:site] + "1" + two[site+1:]
            else:
                one = one[:site] + "1" + one[site+1:]
            new = new[:site] + "1" + new[site+1:]
            # change the "0" that is at the site to a "1"
        else:
            # if the flag is true that means that all of the sites have been mutated,so the allele string stays the same
            one = one
            two = two
            new = new
    t.set_node_property('copy 1 alleles', one)
    t.set_node_property('copy 2 alleles', two)
    t.set_node_property('alleles', new)
    # set the allele string of the cell
    t.set_node_property('new mutations', mutated_sites)
    # setting the new mutations to be used to build the mutation tree
    if t.is_leaf() == False:
        set_alleles(t.get_left_child(), sites, t)
        set_alleles(t.get_right_child(),sites, t)

def create_neutral_phylogenetic(numcells, numsites, percent_clonal):
    """
    :param numclones: The number of clones that will be at the leaves of the tree
    :return: A phylogenetic tree that models branching cancer cell evolution in the form of a Newick string
    """
    sequences = []
    total_length = 0
    changed_sites = []
    for int in range(numcells):
        sequences.append("cell" + str(int + 1))
    trees = []
    num2 = 0
    num3 = 0
    total = 0
    for sequence in sequences:
        num = FullBiTree(sequence)
        #creates a node for each cell that will be leaf
        trees.append(num)
        #adds this node to the trees array
        num2 += 1
    while len(trees) > 1:
        if len(trees) == 1:
            return trees[0]
        num3 += 1
        random1 = random.randint(0, len(trees) - 1)
        random2 = random.randint(0, len(trees) - 1)
        # choosing two random leaves
        if random1 == random2:
            while random2 == random1:
                random2 = random.randint(0, len(trees) - 1)
        lnum = random.random()
        rnum = random.random()
        # the branch lengths of the leaves
        total_length = total_length + lnum + rnum
        left = trees[random1]
        right = trees[random2]
        left.set_node_property('branch length', lnum)
        right.set_node_property('branch length', rnum)
        temp = FullBiTree("Intermediate" + str(num3), left, right)
        #create an intermediate node with these two leaves as the children
        trees.remove(trees[random1])
        if random2 != 0 and random2 > random1:
            trees.remove(trees[random2 - 1])
        else:
            trees.remove(trees[random2])
        # remove the two leaves from being candidates
        trees.append(temp)
        # add the intermidiate node to the list that can be chosen
        num2 += 1
    # l = random.random()
    # total_length += l
    trees[0].set_node_property('branch length', total_length / percent_clonal - total_length)
    normal = FullBiTree("normal")
    normal.set_node_property('alleles', "0" * numsites)
    normal.set_node_property('copy 1 alleles', "0" * numsites)
    normal.set_node_property('copy 2 alleles', "0" * numsites)
    normal.set_node_property('nummutes', 0)
    normal.set_node_property('new mutations', [])
    new = FullBiTree("root", normal, trees[0])
    new.set_node_property('alleles', "0" * numsites)
    new.set_node_property('copy 1 alleles', "0" * numsites)
    new.set_node_property('copy 2 alleles', "0" * numsites)
    new.set_node_property('branch length', 0)
    new.set_node_property('nummutes', 0)
    new.set_node_property('new mutations', [])
    new.set_node_property('clone', new.get_name())
    define_nummutes(new.get_right_child(), total_length / percent_clonal, numsites)
    set_alleles(new.get_right_child(), changed_sites, new)
    # initial conditions for the tree to be biologically accurate and to pass down to children
    if len(changed_sites) != numsites:
        print "REDOING-------------------------"
        new = create_neutral_phylogenetic(numcells, numsites)
    return new

def find_clones(t, alleles = {}, new_mutes = {}, clones = [], clone_leaves = [], copy_one = {}, copy_two = {}):
    """
    :param t: A fullbitree
    :return: The full bitree with all of the cells that have the same clones defined
    """
    left = t.get_left_child().get_node_property('nummutes')
    # if it is 0, then we add it to the clone of the parent
    right = t.get_right_child().get_node_property('nummutes')
    # if it is 0, then we add it to the clone of the parent
    original = t.get_node_property('clone')
    alleles[t.get_node_property('clone')] = t.get_node_property('alleles')
    new_mutes[t.get_node_property('clone')] = t.get_node_property('new mutations')
    copy_one[t.get_node_property('clone')] = t.get_node_property('copy 1 alleles')
    copy_two[t.get_node_property('clone')] = t.get_node_property('copy 2 alleles')
    # the clone of the parent cell
    if left == 0:
        t.get_left_child().set_node_property('clone', original)
        # if there are no mutations in the child node then we set the clone to be the clone of the parent
    elif left != 0:
        t.get_left_child().set_node_property('clone', t.get_left_child().get_name())
        clones.append(t.get_left_child().get_node_property('clone'))
        alleles[t.get_left_child().get_node_property('clone')] = t.get_left_child().get_node_property('alleles')
        new_mutes[t.get_left_child().get_node_property('clone')] = t.get_left_child().get_node_property('new mutations')
        copy_one[t.get_left_child().get_node_property('clone')] = t.get_left_child().get_node_property('copy 1 alleles')
        copy_two[t.get_left_child().get_node_property('clone')] = t.get_left_child().get_node_property('copy 2 alleles')
        # if there are mutations in the child node then create a new clone in the tree
    if right == 0:
        # if there are no mutations in the child node then we set the clone to be the clone of the parent
        t.get_right_child().set_node_property('clone', original)
    elif right != 0:
        # if there are mutations in the child node then create a new clone in the tree
        t.get_right_child().set_node_property('clone', t.get_right_child().get_name())
        clones.append(t.get_right_child().get_node_property('clone'))
        alleles[t.get_right_child().get_node_property('clone')] = t.get_right_child().get_node_property('alleles')
        new_mutes[t.get_right_child().get_node_property('clone')] = t.get_right_child().get_node_property('new mutations')
        copy_one[t.get_right_child().get_node_property('clone')] = t.get_right_child().get_node_property('copy 1 alleles')
        copy_two[t.get_right_child().get_node_property('clone')] = t.get_right_child().get_node_property('copy 2 alleles')
    if t.get_left_child().is_leaf() == True:
        if t.get_left_child().get_node_property('clone') not in clone_leaves:
            #print "FLAG"
            clone_leaves.append(t.get_left_child().get_node_property('clone'))
    if t.get_right_child().is_leaf() == True:
        if t.get_right_child().get_node_property('clone') not in clone_leaves:
            #print "FLAG"
            clone_leaves.append(t.get_right_child().get_node_property('clone'))
    if t.get_left_child().is_leaf() != True:
        find_clones(t.get_left_child(), alleles, new_mutes, clones, clone_leaves, copy_one, copy_two)
    if t.get_right_child().is_leaf() != True:
        find_clones(t.get_right_child(), alleles, new_mutes, clones, clone_leaves, copy_one, copy_two)
    # recursively call it on the children of the cell we are currently at in the tree

def convert_to_mutation(t, root):
    """
    :param t: A phylogenetic tree that the user wants to convert
    :param root: The root for the mutation tree
    :return: A mutation tree that represents the phylogenetic tree passed into this function
    """
    left_mutes = t.get_left_child().get_node_property('new mutations')
    right_mutes = t.get_right_child().get_node_property('new mutations')
    # the new mutations in both of the children
    common_mutes = []
    # will store the mutations that are shared between the left and right
    for mute1 in left_mutes:
        for mute2 in right_mutes:
            if mute1 == mute2:
                # this is never going to happen but just in case?
                common_mutes.append(mute1)
    temp = root
    for mutation in common_mutes:
        # there will never be anything in common mutes, idk why this code is still here
        a = Tree(str(mutation))
        temp.add_children([a])
        temp = a
    new = temp
    for mutation in left_mutes:
        if mutation not in right_mutes:
            # if there are multiple mutations in one branch, this will stack them
            a = Tree(str(mutation))
            temp.add_children([a])
            temp = a
    if t.get_left_child().is_leaf() == False:
        # only recursively call when the children of the node you are at are not leaves
        convert_to_mutation(t.get_left_child(), temp)
    for mutation in right_mutes:
        if mutation not in left_mutes:
            # if there are multiple mutations in one branch this will stack them
            a = Tree(str(mutation))
            new.add_children([a])
            new = a
    if t.get_right_child().is_leaf() == False:
        # only recursively call when the children are not leaves
        convert_to_mutation(t.get_right_child(), new)

def can_mutate(t, poss):
    """
    :param t: a mutation tree
    :param poss: the possible progressions
    :return: a dictionary containing which mutations lead to eachother
    """
    for mutation in t.get_children():
        poss[t.get_name()].append(int(mutation.get_name()))
        can_mutate(mutation, poss)

def find_fitness(mean, stdev, passprob, numsites, t):
    """
    :param numsites: the number of sites that can be changed
    :param br: the birth rate of the initial cell
    :param dr: teh death rate of the initial cell
    :param mean: the mean of the normal distribution that the added fitness will be sampled from
    :param stdev: the standard deviation of the normal distribution that the added fitness will be sampled from
    :param passprob: the probability that a mutation does not add any fitness
    :return: each of the mutations will now be assigned a fitness that will be used to determine overall size at the
     end of the simulation, and also which clones can evolve into which other clones will be determined
    """
    mute_tree = Tree("root")
    convert_to_mutation(t, mute_tree)
    # the mutation tree of teh phylogenetic tree passed into the function
    all = defaultdict(str)
    # the alleles of the clones
    clones = []
    # all of the clones
    new_mutes = defaultdict(list)
    # each clone mapped to its new mutations
    fit = defaultdict(float)
    # each site mapped to its added fitness
    clone_fit = defaultdict(float)
    # the added fitness of each clone
    clone_leaves = []
    copy_one = defaultdict(str)
    copy_two = defaultdict(str)
    find_clones(t, all, new_mutes, clones, clone_leaves, copy_one, copy_two)
    poss = defaultdict(list)
    # the mutations mapped to all of the mutations that they can evolve into
    can_mutate(mute_tree, poss)
    #print poss
    clone_poss = defaultdict(list)
    for site in range(numsites):
        ran = random.random()
        if ran > passprob:
            added = np.random.normal(mean, stdev)
            fit[site] = added
    used_clones = []
    for clone in clones:
        for num in range(len(all[clone])):
            if all[clone][num] == '1':
                clone_fit[clone] += fit[num]
                if num in new_mutes[clone]:
                    for clone2 in clones:
                        if all[clone2][num] == '1' and clone2 != clone:
                            clone_poss[clone].append(clone2)
                            used_clones.append(clone2)
    for clone in clones:
        if clone not in used_clones:
            clone_poss["root"].append(clone)
    clones.append("root")
    return clone_poss, clone_fit, clones

def find_clone_sizes(clone_poss, clone_fit, clones, br, dr, mr, total_count):
    clone_dr = defaultdict(float)
    clone_cells = defaultdict(int)
    for clone in clones:
        # sets the initial birthrate for all clones
        clone_fit[clone] += br
        clone_dr[clone] += dr
    birthing = defaultdict(None)
    for clone in clones:
        # a dictionary of booleans for each clone to see if it is reproducing or not
        birthing[clone] = False
    birthing["root"] = True
    # sets the original clone to be true
    cell_count = 1
    clone_cells["root"] += 1
    # adds a cell for teh root clone
    while cell_count < total_count:
        #does this process until the maximum number of cells is reached
        for clone in clones:
            if birthing[clone] == True:
                # goes through all of the clones taht are reproducing
                count = int(clone_cells[clone])
                for num1 in range(count):
                    # for each cell that is alive
                    for num in range(int(round(clone_fit[clone]))):
                        rand = random.random()
                        if rand < mr and clone_poss[clone] != []:
                            # randomly mutates
                            new_clone = clone_poss[clone][random.randint(0,len(clone_poss[clone])-1)]
                            clone_poss[clone].remove(new_clone)
                            birthing[new_clone] = True
                            clone_cells[new_clone] += 1
                            cell_count += 1
                        else:
                            # if there is no mutation just reproduce into the same clone
                            clone_cells[clone] += 1
                            cell_count += 1
        for clone in clones:
            if birthing[clone] == True and clone_cells[clone] > 1:
                # if the clone is alive kill off the death rate number of cells
                clone_cells[clone] -= clone_dr[clone]
                cell_count -= clone_dr[clone]
    return clone_cells

def find_percentages(clone_cells):
    """
    :param clone_cells: the number of cells in each clone
    :return: the percentages of the total cells in the tumor that each clone contains
    """
    total_count = 0.0
    percentages = defaultdict(float)
    for clone in clone_cells:
        total_count += clone_cells[clone]
    for clone in clone_cells:
        percentages[clone] = float(clone_cells[clone]) / total_count
    return percentages

def mutate_genome(genome, t):
    """
    :param genome: the actual string of the genome "A", "C", "T", "G", etc
    :param leaves: the leaves of the clonal tree that we will be transforming
    :return: the genome for each of the clones at the leaves of the clonal tree
    """
    used_sites = []
    all = defaultdict(str)
    copy_one = defaultdict(str)
    copy_two = defaultdict(str)
    # the alleles of the clones
    clones = []
    # all of the clones
    new_mutes = defaultdict(list)
    # each clone mapped to its new mutations
    clone_leaves = []
    upper = ["A", "C", "T", "G"]
    lower = ["a" , "c", "t", "g"]
    find_clones(t, all, new_mutes, clones, clone_leaves, copy_one, copy_two)
    #finds the clones in the tree
    genome_change = defaultdict(str)
    # will hold the mutations
    genome_copies = defaultdict(list)
    sites = defaultdict(int)
    for num in range(len(copy_one[clones[0]])):
        # goes through the first copy
         rand = random.randint(0,len(genome) - 1)
        #randomizes a spot in the genome
         while rand in used_sites:
             #changes the spot until it has not been changed it
             rand = random.randint(0, len(genome) - 1)
         used_sites.append(rand)
         sites[num] = rand
    for num in range(len(copy_one[clones[0]])):
        # stores what the sight will be changed to
        if genome[sites[num]] in lower:
            temp = copy.copy(lower)
            temp.remove(genome[sites[num]])
            rand = random.randint(0,2)
            genome_change[num] = temp[rand]
        elif genome[sites[num]] in upper:
            temp = copy.copy(upper)
            temp.remove(genome[sites[num]])
            rand = random.randint(0,2)
            genome_change[num] = temp[rand]
    change = {}
    for site in sites:
        #stores the actual changed site
        change[sites[site]] = (genome[sites[site]], genome_change[site])
    clone_leaves.remove("root")
    for clone in clone_leaves:
        # actually changes the site
        one = copy_one[clone]
        two = copy_two[clone]
        first = copy.deepcopy(genome)
        second = copy.deepcopy(genome)
        for num in range(len(one)):
            if one[num] == "1":
                first = first[:sites[num]] + genome_change[num] + first[sites[num]+1:]
        for num in range(len(two)):
            if two[num] == "1":
                second = second[:sites[num]] + genome_change[num] + second[sites[num] + 1:]
        copies = [first, second]
        genome_copies[clone] = copies
    return genome_copies, sites, change

def sample_cells(percentages, clones, num_cells):
    nump = defaultdict(int)
    total = 0
    for clone in clones:
        nump[clone] = int(percentages[clone] * float(num_cells))
        total += int(percentages[clone] * float(num_cells))
    while total < num_cells:
        nump[random.choice(clones.keys())] += 1
        total += 1
    while total > num_cells:
        nump[random.choice(clones.keys())] -= 1
        total -= 1
    sampled = defaultdict(list)
    for clone in clones:
        for num in range(nump[clone]):
            sampled[clone+"_"+str(num)] = clones[clone]
    return sampled

def amplify_genome(genome_copies, fp_rate, num_rounds, ado_rate, num_sections):
    """
    :param cell_dict: the dictionary that maps the leaves of the tree to its two cell copies
    :param drop_rate: the probability that a spot in the genome gets dropped
    :param fp_rate: the probability that a spot in the genome gets changed incorrectly
    :return: fasta files of the copied parts of the genome
    """
    upper = ["A", "C", "T", "G"]
    lower = ["a", "c", "t", "g"]
    # p: Dictionary of probabilities
    # c: Dictionary of counts
    p = defaultdict(lambda: defaultdict(lambda: defaultdict()))
    c = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 1)))
    sections = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: [])))
    # muts = defaultdict(list)
    # NUM_UNITS = 10
    for clone in genome_copies:
        # Make sure the two copies are of the same length
        min_len = min(len(genome_copies[clone][0]), len(genome_copies[clone][1]))
        cur_pos = 0
        unit_length = min_len / num_sections
        for j in range(1, num_sections):
            # Use Normal distribution to chunk up the genomes
            end_pos = int(random.gauss(j*unit_length, 20))
            sections[clone][j][0].append(genome_copies[clone][0][cur_pos:end_pos])
            sections[clone][j][1].append(genome_copies[clone][1][cur_pos:end_pos])
            cur_pos = end_pos
            rand3 = random.random()
            # Set the initial probabilities for two copies
            if rand3 < ado_rate:
                rand4 = random.random()
                if rand4 < 0.5:
                    p[clone][j][0] = 0.9
                    p[clone][j][1] = 0.1
                else:
                    p[clone][j][1] = 0.9
                    p[clone][j][0] = 0.1
            else:
                p[clone][j][0] = 0.5
                p[clone][j][1] = 0.5
            for num in range(num_rounds):
                rand4 = random.random()
                # Update copy counts and random selection of genome copy.
                if rand4 < p[clone][j][0]:
                    flag = 0
                    c[clone][j][0] += 1
                    amp = random.choice(sections[clone][j][0])
                else:
                    flag = 1
                    c[clone][j][1] += 1
                    amp = random.choice(sections[clone][j][1])
                # Update probability if ADO occurs
                if rand3 < 0.2:
                    p_temp_0 = p[clone][j][0]
                    p_temp_1 = p[clone][j][1]
                    p[clone][j][0] = c[clone][j][0] * p_temp_0 / (c[clone][j][0] * p_temp_0 + c[clone][j][1] * p_temp_1)
                    p[clone][j][1] = c[clone][j][1] * p_temp_1 / (c[clone][j][0] * p_temp_0 + c[clone][j][1] * p_temp_1)
                # Make a copy of the chosen section
                place = random.randint(0, len(amp)/10)
                temp = amp[place:]
                # Introducing FP errors
                for num2 in range(len(temp)):
                    rand = random.random()
                    if rand < fp_rate:
                        if temp[num2] in upper and temp[num2] != "N":
                            temp2 = deepcopy(upper)
                            temp2.remove(temp[num2])
                            mut = random.choice(temp2)
                            temp = temp[:num2] + mut + temp[num2 + 1:]
                            # muts[clone].append((num2, temp[num2], mut))
                        elif temp[num2] in lower and temp[num2] != "n":
                            temp2 = deepcopy(lower)
                            temp2.remove(temp[num2])
                            mut = random.choice(temp2)
                            temp = temp[:num2] + mut + temp[num2 + 1:]
                            # muts[clone].append((num2, temp[num2], mut))
                        # print "temp num2 after", temp[num2]
                sections[clone][j][flag].append(temp)
                # sections[clone][j][flag][0] += temp
    return sections
