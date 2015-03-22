##
## Step 2b - exclusion of redundant reads
##

file 'tmp/step2b' => 'tmp/step2a' do |t|
  outfp = open("| pigz -c > #{t.name}", 'w')
  tracefp = open("| gsort --parallel=#{PROCS} -S #{3*PROCS}G -k 1,1 | pigz -c > #{t.name}.trace", 'w')
  infp = open("| unpigz -c #{t.source}")
  line = infp.gets
  prebc, preacc, preqv, preseq = line.rstrip.split(/\t/)
  tracefp.puts [prebc, preacc, preacc].join("\t")
  while line = infp.gets
    bc, acc, qv, seq = line.rstrip.split(/\t/)
    unless preseq == seq
      outfp.puts [prebc, preacc, preqv, preseq].join("\t")
      prebc = bc
      preacc = acc
      preqv = qv
      preseq = seq
    end
    tracefp.puts [acc, preacc].join("\t")
  end
  infp.close
  tracefp.close
  outfp.close
end

task :clean_step2b do
  rm_rf "tmp/step2b"
  rm_rf "tmp/step2b.trace"
end
