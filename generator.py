import gzip
import os
from trees import *
from generator_helper import *
from genome_helper import *


def main(num_clones, num_cells, num_mutations, fp_rate, ado_rate, ref_file, num_sections):
    # Take out spaces in the genome
    # chr20 = open("chr20.fa","r")
    ref_genome = open(ref_file, "r")
    new_genome = ref_genome.read()
    new_genome_2 = new_genome.replace('\n', '')
    new_genome = new_genome_2[66341:166341]
    # Simulate single cell genomes
    neut = create_neutral_phylogenetic(num_clones, num_mutations)
    a = mutate_genome(new_genome, neut)
    # print "genome copies", a[0]
    # print "sites", a[1]
    # print "change", a[2]
    f = find_fitness(4, .15, .02, num_mutations, neut)
    s = find_clone_sizes(f[0], f[1], f[2], 2.0, 2.0, .01, 1000000)
    p = find_percentages(s)
    sampled = sample_cells(p, a[0], num_cells)
    # Build the matrix that identifies mutations
    matrix = [["pos", "ref","alt"]]
    for pos, bases in a[2].items():
        pos3 = pos + 66336
        matrix.append([pos3, bases])
    for cellname, genome in sampled.items():
        matrix[0].append(cellname)
        row = 1
        for pos2, bases2 in a[2].items():
            if bases2[0] == " ":
                matrix[row].append(2)
            # Append 0 if it is the ref genome
            elif genome[0][pos2] == bases2[0]:
                matrix[row].append(0)
            # Append 1 if it is the alt genome
            elif genome[0][pos2] == bases2[1]:
                matrix[row].append(1)
            else:
                matrix[row].append(3)
            row += 1
    # Write the matrix to a txt file
    f = open("mutation.txt", "w")
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            f.write(str(matrix[i][j]) + " ")
        f.write("\n")
    f.close()
    amps = amplify_genome(sampled, fp_rate, 10, ado_rate, num_sections)
    err_file = "errorModel/ill100v4_s.gzip"
    # err_file = "D:/2018 summer/Monovar/GemSIM_v1.6/models/ill100v4_s.gzip"
    # For each cell, write all genome reads to a single file
    ls_fq = open("list_fastq.txt","w")
    for cell in amps:
        ls_fq.write(str(cell)+"\n")
        for unit in amps[cell]:
            for amp in amps[cell][unit]:
                for idx in range(len(amps[cell][unit][amp])):
                    fasta = amps[cell][unit][amp][idx]
                    # print len(fasta)
                    # print "fasta",fasta
                    GenReads(fasta=fasta, length=100, qual=33, models=err_file, out=str(cell))
                    # out= str(cell) + "s" + str(unit) + "ADO" + str(amp) + "cp" + str(idx))
                    # break
                # break
            # break
        # break


main(num_clones=10, num_cells=100, num_mutations=50, fp_rate=10**(-4),
     ado_rate=0.2, ref_file="chr20.fa", num_sections=10)
