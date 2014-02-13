require 'parallel'

load 'rake/configuration.rake'
load 'rake/removePhyX.rake'
load 'rake/demultiplex.rake'

task :default => :demultiplex

##

# def read_barcodes
#   wells = Array.new
#   open('src/barcodes.txt').each do |line|
#     wells.push(line.rstrip.split[0])
#   end
#   wells
# end

# WELLS = read_barcodes

# ##

# task :default => [:alignment]

# ##

# task :alignment => ['out/stat/alignment.txt', 'out/stat/spike.txt.gz']

# BAMS = Array.new
# WELLS.each { |well|
#   bam = "out/ali/#{LIBID}.#{well}/accepted_hits.bam"
#   file bam => ["out/seq/#{LIBID}.#{well}.fq.xz", 'src/knownGene.annotated.gtf', 'src/ebwt/ref.1.ebwt'] do |t|
#     align(t)
#   end
#   BAMS.push(bam)
# }

# def align(t)
#   outdir = t.name.sub('/accepted_hits.bam', '')
#   logfile = t.name.sub('out/ali/', 'log/align.').sub('/accepted_hits.bam', '.log')
#   gtf = Dir.exist?('tmp/knownGene.1.ebwt') ? '' : "-G #{t.prerequisites[1]}"
#   sh "mkdir -p #{outdir}"
#   sh <<"PROCESS"
# xzcat #{t.prerequisites[0]} | \\
# tophat #{gtf} --transcriptome-index tmp/knownGene --num-threads #{PROCS} \\
#        --library-type fr-secondstrand --min-anchor 5 --coverage-search \\
#        --bowtie1 --output-dir #{outdir} \\
#        #{t.prerequisites[2].sub('.1.ebwt', '')} /dev/stdin > #{logfile} 2>&1
# PROCESS
# end

# file 'out/stat/spike.txt.gz' => BAMS do |t|
#   fp = open("| gzip -c --best > #{t.name}", 'w')
#   t.prerequisites.each { |bam|
#     head = [LIBID, /#{LIBID}\.([^\/]+)/.match(bam).to_a[1], ''].join("\t")
#     open("| samtools view #{bam} | cut -f 1,3,4,12- | grep '\tRNA_SPIKE_' | grep XS:A:+ | cut -f 2,3").each { |line| fp.puts head + line }
#   }
#   fp.close
# end

# file 'out/stat/alignment.txt' => BAMS do |t|
#   well2cnts = Hash.new
#   Parallel.each(t.prerequisites, :in_threads => PROCS) { |mbam|
#     tmp = Array.new
#     fp = open(mbam.sub('accepted_hits.bam', 'align_summary.txt'))
#     fp.gets
#     tmp.push(/Input\:\s+(\d+)/.match(fp.gets).to_a[1])
#     mapped = /Mapped\:\s+(\d+)/.match(fp.gets).to_a[1]
#     multi = /these\:\s+(\d+)/.match(fp.gets).to_a[1]
#     tmp.push(mapped.to_i-multi.to_i)
#     tmp.push(mapped)
#     fp.close
#     open("| samtools view #{mbam} | cut -f 1,3,12- | grep '\tRNA_SPIKE_' | grep XS:A:+ | cut -f 1 | sort -u | wc -l").each { |line|
#       tmp.push(line.rstrip) if line != "\n"
#     }
#     open("| samtools view #{mbam} | cut -f 1,3 | grep -E '\t(U13369|RIBOSOMAL)' | cut -f 1 | sort -u | wc -l").each { |line|
#       tmp.push(line.rstrip) if line != "\n"
#     }
#     well = /#{LIBID}\.([^\/]+)/.match(mbam).to_a[1]
#     well2cnts[well] = [LIBID, well] + tmp
#   }
#   # open("| samtools view #{mbam} | cut -f 1 | wc -l").each { |line|
#   #   well2cnts[well].push(line.rstrip) if line != "\n"
#   # }
#   fp = open(t.name, 'w')
#   fp.puts ['LIB', 'WELL', 'TOTAL', 'MAPPED.UNIQUE', 'MAPPED.ALL', 'SPIKE', 'RIBOSOMAL'].join("\t")
#   well2cnts.keys.sort.each { |well| fp.puts well2cnts[well].join("\t") }
#   fp.close
# end
