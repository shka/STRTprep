tmp = Array.new
LIBIDS.each { |libid|
  taskid = "removePhyX_#{libid}"
  tmp.push(taskid)
  task taskid => "out/ali/#{libid}.phyX.bam"
}

task 'removePhyX' => tmp

rule '.phyX.bam' => proc { |target|
  libid = target.pathmap("%n").pathmap("%n")
  tmp = [File.expand_path(CONF['LIBS'][libid]['PHYX']+".1.ebwt")]
  CONF['LIBS'][libid]['FASTQS'].each { |fastq|
    tmp.push(File.expand_path(fastq))
  }
  tmp
} do |t|
  sh 'mkdir -p out/ali log'
  libid = t.name.pathmap("%n").pathmap("%n")
  sh "(gunzip -c #{t.prerequisites[1..-1].join(' ')} | bowtie -S -p #{PROCS} -t #{t.prerequisites[0].sub('.1.ebwt', '')} - | samtools view -@ #{PROCS} -S -b -o #{t.name} -) 2> log/#{libid}.removePhyX.log "
end
