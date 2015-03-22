##
## Step 3a - Decollapse sequences & alignments
##

rule /.step3a/ => ['tmp/step2a', 'tmp/step2b.trace'] do |t|
  repacc2accpath = mymkfifo('step3a-')
  pid1 = spawn "unpigz -c #{t.sources[1]} > #{repacc2accpath}"

  libwellid = t.name.pathmap("%n")
  bcaccpath = mymkfifo('step3a-')
  pid2 = spawn <<EOF
unpigz -c #{t.source} \
| grep '^#{libwellid}\t'\
| gcut -f 2-4 \
| gsort --parallel=#{PROCS} -S #{1500*PROCS}M -k 1,1 > #{bcaccpath}
EOF
# | parallel --pipe -L 10000 --round-robin grep '^#{libwellid}\t'\

  sh <<EOF
gjoin -t '\t' -j 1 -o 1.2,2.1,2.2,2.3 #{repacc2accpath} #{bcaccpath}\
| gsort --parallel=#{PROCS} -S #{1500*PROCS}M -k 1,1\
| pigz -c > #{t.name}
EOF

  Process.waitpid(pid1)
  Process.waitpid(pid2)
end

#

def step3a_bam_sources(path)
  return ['tmp/step2c/accepted_hits.header',
          'tmp/step2c/accepted_hits.samSortedByAcc',
          path.sub(/^out\/bam/, 'tmp').sub(/\.bam$/, '.step3a')]
end

rule /^out\/bam\/[^\/]+\.bam$/ => [->(path){ step3a_bam_sources(path) }] do |t|
  mkdir_p t.name.pathmap('%d')
  
  outfifo = mymkfifo('step3a-bam-')
  pid = spawn "samtools view -@ #{PROCS} -b #{outfifo} > #{t.name}"

  tmpfifo1 = mymkfifo('step3a-bam-')
  pid1 = spawn "gunzip -c #{t.sources[1]} > #{tmpfifo1}"

  tmpfifo2 = mymkfifo('step3a-bam-')
  pid2 = spawn "gunzip -c #{t.sources[2]} > #{tmpfifo2}"

  outfp = open(outfifo, 'w')
  infp = open(t.source)
  while line = infp.gets
    outfp.puts line
  end
  infp.close
  infp = open("| gjoin -t '\t' -j 1 #{tmpfifo1} #{tmpfifo2} | gsort --parallel=#{PROCS} -S #{3*PROCS} -k 3,3 -k 4,4n")
  while line = infp.gets
    cols = line.rstrip.split(/\t/)
    cols[0] = cols[-3]
    unless cols[11..-4].include?('XS:A:-')
      cols[9] = cols[-1]
      cols[10] = cols[-2]
    else
      cols[9] = cols[-1].reverse.tr('ACGT', 'TGCA')
      cols[10] = cols[-2].reverse
    end
    outfp.puts cols[0..-4].join("\t")
  end
  infp.close
  outfp.close

  Process.waitpid(pid)
  Process.waitpid(pid1)
  Process.waitpid(pid2)
end

#

task :clean_step3a do
  LIBIDS.each do |libid|
    rm_rf "#{libid}.*.step3a"
  end
end
