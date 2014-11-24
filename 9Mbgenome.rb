#encoding: utf-8
require_relative 'lib/model_genome'
require_relative 'lib/write_it'
name = ARGV[0]
require 'PMeth'

# make the directory to put data files into
Dir.mkdir(File.join(Dir.home, "/Correlations/genomes/#{name}"))

# Create the lists of homozygous and heterozygous SNPs
hm_r = 'hm <- rnorm(9000, 4500000, 1000)' # Causative SNP at/near 10000
ht_r = 'ht <- runif(9000, 1, 9000000)'   # Genome length of 10000
hm, ht = ModelGenome::get_snps(hm_r, ht_r)
snp_pos = [hm, ht].flatten

puts "There are #{hm.length} homozygous SNPs"
puts "There are #{ht.length} heterozygous SNPs"
puts "Is there a SNP at the centre of the distribution? -- #{snp_pos.include?(500000)}"

arabidopsis_c4 = ModelGenome::fasta_to_char_array("TAIR10_chr4.fasta")
puts "Creating the genome..."
small_genome = arabidopsis_c4[-9000000..-1] # Genome length of 100 kb
contig_size = 9000 # 100-200 bp
puts "...and generating the fragments"
frags = ModelGenome::get_frags(small_genome, contig_size)

puts "Small genome     length: #{small_genome.length} bases"
puts "You have created #{frags.length} fragments of sizes #{contig_size}-#{contig_size*2}"


# Get the positions of the SNPs on fragments
pos_on_frags, snp_pos_all = ModelGenome::pos_each_frag(snp_pos, frags)

fastaformat_array = ModelGenome::fasta_array(frags)

#Shuffle the original order of fragments by using the adjacent_swap method. 
start_pop = []
start_pop << fastaformat_array
10.times do
	adj_pop = []
	start_pop.each do |perm|
		70.times do
			perm = PMeth.adjacent_swap(perm) 
		end
		adj_pop << perm # need this population to be the next starting population
	end
	start_pop = adj_pop
end

fastaformat_array_shuf = start_pop
vcf = ModelGenome::vcf_array(frags, pos_on_frags, snp_pos_all, hm, ht)

WriteIt::write_data("genomes/#{name}/frags.fasta", fastaformat_array)
WriteIt::write_data("genomes/#{name}/snps.vcf", vcf)
WriteIt::write_data("genomes/#{name}/frags_shuffled.fasta", fastaformat_array_shuf)
WriteIt::write_txt("genomes/#{name}/info", [hm_r, ht_r, "Contig size = #{contig_size}"])
WriteIt::write_txt("genomes/#{name}/hm_snps", hm)
WriteIt::write_txt("genomes/#{name}/ht_snps", ht)