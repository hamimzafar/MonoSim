import subprocess

acclist_file = open("list_fastq.txt", "r")
acclist = acclist_file.read().strip().split()

print("Acclist: ", acclist)

# Download fasta files
try:
    subprocess.run(["mkdir", "fq"])
except subprocess.CalledProcessError:
    print("Directory already created.")

for seq in acclist:
    # Align
    with open("fq/" + seq + ".sam", "w") as outfile:
        subprocess.call(["bwa", "mem", "chr20.fa", "fq/" + seq + ".fastq"], stdout=outfile)
    print("Saved as " + seq + ".sam")

    # Sort samfile
    subprocess.check_output(
        ["samtools", "sort", "-o", "fq/" + seq + "_sorted.bam", "fq/" + seq + ".sam"])

    # Cleanup
    subprocess.run(["rm", "fq/" + seq + ".fastq"])
    subprocess.run(["rm", "fq/" + seq + ".sam"])