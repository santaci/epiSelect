import pysam
import argparse
import os
import subprocess

parser = argparse.ArgumentParser(description='Make VCFs')
parser.add_argument('-f','--file', help='msprime VCF', required=True)
parser.add_argument('-c','--chrom', help='chromosome', required=True, type=int)
parser.add_argument('-g1','--gen1', help='generation 1', required=True, type=int)
parser.add_argument('-g2','--gen2', help='generation 2', required=False, type=int)
parser.add_argument('-bp','--breakpoints', help='breakpoints', required=True, nargs="*")
args = parser.parse_args()
print(args.breakpoints)
# Define input and output VCF file paths
input_vcf = args.file + ".vcf.gz"  # Replace with your actual VCF.gz file path

if not os.path.exists(input_vcf+".tbi"):
    pysam.tabix_index(input_vcf,preset='vcf')

output_vcf_prefix = args.file
g1 = str(args.gen1)
if args.gen2 is None:
    g2 = 0
else:
    g2 = str(args.gen2)
# Define nucleotide assignment for REF and ALT
allele_dict = {0: 'A', 1: 'T'}  # Example assignment: 0 -> 'A', 1 -> 'T'

# Open the input VCF
vcf_in = pysam.VariantFile(input_vcf, "r")

# Add the 'AA' INFO field to the header if it doesn't exist
if 'AA' not in vcf_in.header.info:
    vcf_in.header.info.add("AA", 1, "String", "Ancestral Allele")

for contig, metadata in vcf_in.header.contigs.items():
    total_length = metadata.length
    print(f"Original Contig: {contig}, Length: {metadata.length}")
    contig = int(args.chrom)

# Define breakpoints and new chromosome names
breakpoints = args.breakpoints[0].split(' ')
breakpoints = [int(bp) for bp in breakpoints]
new_chrom_names = [f"chr{i+contig}" for i in range(len(breakpoints) + 1)]

# Calculate lengths of new chromosomes based on breakpoints
chrom_lengths = []
chrom_lengths.append(breakpoints[0])  # First segment length
for i in range(1, len(breakpoints)):
    chrom_lengths.append(breakpoints[i] - breakpoints[i - 1])
chrom_lengths.append(total_length - breakpoints[-1])  # Last segment length goes to the end of the chromosome
print(chrom_lengths)

# Function to update the header with new chromosome names and lengths
def update_header(vcf_header, chrom_names, chrom_lengths):
     # Add new contigs (chromosomes) with their lengths
    for chrom, length in zip(chrom_names, chrom_lengths):
        if chrom not in list((vcf_header.contigs)):
            vcf_header.contigs.add(chrom, length=length)

    # Create a copy of the existing header
    new_header = vcf_header.copy()
    return new_header

# Create output VCF files for each new chromosome segment
if args.file.__contains__("gwas"):
    print("GWAS")
    combined_vcf = pysam.VariantFile(f"{output_vcf_prefix}_{g1}_{g2}.vcf.gz", "w", header=update_header(vcf_in.header, new_chrom_names, chrom_lengths))
else:
    print("iHS")
    if (g1 != g2) and (g2 == 0):
         vcf_outs = {
        name: pysam.VariantFile(f"{output_vcf_prefix}_{name}_d.vcf.gz", "w", header=update_header(vcf_in.header, [name], [chrom_lengths[i]]),) for i, name in enumerate(new_chrom_names)
        }
    else:
        vcf_outs = {
        name: pysam.VariantFile(f"{output_vcf_prefix}_{name}.vcf.gz", "w", header=update_header(vcf_in.header, [name], [chrom_lengths[i]]),) for i, name in enumerate(new_chrom_names)
        }
print(new_chrom_names)



# Function to determine the new chromosome and adjusted position based on breakpoints
def get_new_chrom_and_pos(pos, breakpoints):
    for i, bp in enumerate(breakpoints):
        if pos < bp:
            new_pos = pos if i == 0 else pos - breakpoints[i-1]
            return new_chrom_names[i], new_pos
    return new_chrom_names[-1], pos - breakpoints[-1]

# Function to update REF, ALT, and AA based on allele_dict
def update_alleles(record):
    ref_allele = allele_dict.get(str(record.ref), 'A')  # Default to 'N' if ref is missing
    alt_allele = allele_dict.get(str(record.alts), 'G')  # Default to 'N' if alt is missing
    # Update REF and ALT fields
    record.ref = ref_allele
    record.alts = (alt_allele,)
    # Update AA field to match the REF allele (if AA exists)
    record.info.__setitem__('AA',ref_allele)
    return record


# Process each record in the VCF
for record in vcf_in.fetch():
    new_rec = record.copy()
    # Determine the new chromosome and adjusted position
    new_chrom, new_pos = get_new_chrom_and_pos(record.pos, breakpoints)
    # Update the record with new chromosome, position, and alleles
    new_rec.contig=new_chrom
    new_rec.pos=new_pos
    updated_record = update_alleles(new_rec)
    # Write the updated record to the corresponding output VCF file
    if args.file.__contains__("gwas"):
        combined_vcf.write(updated_record)
    else:
        vcf_outs[new_chrom].write(updated_record)
        #os.remove(args.file + ".inds")

# Close all VCF files
vcf_in.close()
if args.file.__contains__("gwas"):
    combined_vcf.close()
    dirname=os.path.dirname(os.path.abspath(__file__)) 
    #+ "/scripts/"
    gt_vcf_gwas=dirname + "/gt_vcf_gwas.sh"
    process = subprocess.Popen(['bash', gt_vcf_gwas, args.file, "0", "0", g1, g2])
    process.wait()
else:
    for vcf_out in vcf_outs.values():
        vcf_out.close()
os.remove(input_vcf)
os.remove(input_vcf + ".tbi")


print("Modified VCFs have been saved with the prefix:", output_vcf_prefix)
