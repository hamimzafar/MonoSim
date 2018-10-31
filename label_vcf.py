mutation = open("mutation.txt", "r")
lines = mutation.readlines()
first = True
set_mut = set()
for line in lines:
    #print line
    if first:
        first = False
        continue
    else:
        pos = ""
        for cha in line:
            if not cha == "(":
                pos += cha
            else:
                set_mut.add(int(pos))
monovar_result = open("monovar_sim_new.vcf", "r")
new_result = open("monovar_fps.txt", "w")
new_result.write("# line  genome pos  mutation type \n")
res_lines = monovar_result.readlines()
count = 1
num_mut = 0
for line2 in res_lines:
    #print line2
    if count <= 20:
        count += 1
        continue
    else:
        count2 = 1
        #print line2
        pos2 = ""
        for cha2 in line2:
            if count2 <= 6:
                count2 += 1
                continue
            else:
                if not cha2 == "\t":
                    pos2 += cha2
                else:
                    if int(pos2) in set_mut:
                        num_mut += 1
                        new_result.write(str(count) + "       ")
                        new_result.write(pos2 + "      ")
                        new_result.write(" True Mutation")
                        new_result.write("\n")
                    else:
                        new_result.write(str(count) + "      ")
                        new_result.write(pos2 + "      ")
                        new_result.write(" False Positive")
                        new_result.write("\n")
                    count += 1
                    break
print num_mut