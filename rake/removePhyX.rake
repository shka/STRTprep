####
#
# dummy paths
#
tmp = Array.new
LIBIDS.each { |libid| tmp.push("out/ali/#{libid}.phyX.bam") }
task :removePhyX => tmp

####
#
# file (out/ali/#{libid}).phyX.bam => [#{libid.PHYX}.1.ebwt]+#{libid.FASTQS}
#
rule '.phyX.bam' => proc { |target|
  conf = CONF[target.pathmap("%n").pathmap("%n")]
  tmp = [File.expand_path(conf['PHYX']+".1.ebwt")]
  conf['FASTQS'].each { |fastq| tmp.push(File.expand_path(fastq)) }
  tmp
} do |t|
  sh 'mkdir -p out/ali out/stat'
  libid = t.name.pathmap("%n").pathmap("%n")
  sh "(gunzip -c #{t.prerequisites[1..-1].join(' ')} | bowtie -S -p #{PROCS} -t #{t.prerequisites[0].sub('.1.ebwt', '')} - | samtools view -@ #{PROCS} -S -b -o #{t.name} -) 2> out/stat/#{libid}.removePhyX.log"
end

####
#
# cleaning
#
task 'clean_removePhyX' do
  LIBIDS.each { |libid|
    sh "rm -rf out/ali/#{libid}.phyX.bam out/stat/#{libid}.removePhyX.log"
  }
end
