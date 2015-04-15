##
## Step 2c - alignment
##

file 'tmp/step2c/reads.fq.gz' => 'tmp/step2b' do |t|
  mkdir_p t.name.pathmap('%d')
  sh <<EOF
unpigz -c #{t.source}\
 | gawk 'BEGIN { OFS="\\n" }; { print "@"$2,$4,"+",$3 }'\
 | pigz -c > #{t.name} 2> #{t.name}.log
EOF
end

file 'tmp/step2c/accepted_hits.bam' =>
     ['tmp/step2c/reads.fq.gz',
      File.expand_path(CONF[LIBIDS[0]]['GENOMESPIKERIBO'])+'.1.ebwt',
      File.expand_path(CONF[LIBIDS[0]]['TRANSCRIPT'])+'.1.ebwt'] do |t|
  sh <<EOF
(tophat\
   --transcriptome-index #{t.sources[2].pathmap('%X').pathmap('%X')}\
   --library-type fr-secondstrand --min-anchor 5 --coverage-search\
   --output-dir #{t.name.pathmap('%d')} --num-threads #{PROCS} --bowtie1\
   --segment-length #{(CONF[LIBIDS[0]]['CDNA'].to_f/2.0).floor}\
   #{t.sources[1].pathmap('%X').pathmap('%X')} #{t.source})\
   > #{t.name}.log 2>&1
EOF
end

rule '.samUniqSortedByAcc' => '.bam' do |t|
  sh <<EOF
samtools view #{t.source}\
| grep -E 'NH:i:1(\\s|$)'\
| gsort --parallel=#{PROCS} -S #{100/(PROCS+1)}% -k 1,1 \
| pigz -c > #{t.name} 2> #{t.name}.log
EOF
end

rule '.header' => '.bam' do |t|
  sh "samtools view -H #{t.source} > #{t.name} 2> #{t.name}.log"
end
