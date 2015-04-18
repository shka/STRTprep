##
## Step 2b - exclusion of redundant reads
##

file 'tmp/step2b' => 'tmp/step2a' do |t|
  pids = Array.new
  
  outfifopath = mymkfifo('step2b-output-')
  pids.push(spawn "cat #{outfifopath} | pigz -c > #{t.name}")

  tracefifopath = mymkfifo('step2b-trace-')
  pids.push(spawn "gsort --parallel=#{PROCS} -S 75% -t '\t' -k 1,1 #{tracefifopath} | pigz -c > #{t.name}.trace")

  outfp = open(outfifopath, 'w')
  tracefp = open(tracefifopath, 'w')
  infp = open("| unpigz -c #{t.source}")
  line = infp.gets
  prebc, preacc, preqv, preseq = line.rstrip.split(/\t/)
  pre = "#{prebc}\t#{preacc}\t#{preqv}\t#{preseq}\n"
  buf = "#{preacc}\t#{preacc}\n"
  while line = infp.gets
    bc, acc, qv, seq = line.rstrip.split(/\t/)
    unless preseq == seq
      if pre.length > 65536
        outfp.write pre
        pre = ''
      end
      pre << "#{bc}\t#{acc}\t#{qv}\t#{seq}\n"
      preacc = acc
      preseq = seq
    end
    if buf.length > 65536
      tracefp.write buf
      buf = ''
    end
    buf << "#{acc}\t#{preacc}\n"
  end
  infp.close
  tracefp.write buf if buf != ''
  tracefp.close
  outfp.write pre if pre != ''
  outfp.close

  pids.each { |pid| Process.waitpid(pid) }
end

file 'tmp/step2b.trace' => 'tmp/step2b'
