##
## Step 2b - exclusion of redundant reads
##

step2b_sources = Array.new
LIBIDS.each do |libid|
  step2b_sources.push("tmp/#{libid}.step2a")
end

file 'tmp/step2b' => step2b_sources do |t|
  pid = fork do
    outfp = open("| pigz -c > #{t.name}", 'w')
    tracefp = open("| gsort --parallel=#{PROCS} -S #{50/(PROCS+1)}% -t '\t' -k 1,1 | pigz -c > #{t.name}.trace", 'w')
    infp = open("| unpigz -c #{t.sources.join(' ')} | gsort --parallel=#{PROCS} -S #{50/(PROCS+1)}% -t '\t' -k 4,4 -k 3,3r")
    preseq = ''
    preacc = ''
    out = ''
    trace = ''
    while line = infp.gets
      bc, acc, qv, seq = line.rstrip.split(/\t/)
      if preseq != seq
        out << "#{bc}\t#{acc}\t#{qv}\t#{seq}\n"
        while 0 < out.bytesize
          begin
            written = outfp.write_nonblock(out)
            out = out.byteslice(written..-1)
          rescue IO::WaitWritable
            if 1048576 < out.bytesize
              IO.select(nil, [outfp])
              retry
            end
          end
        end
        preacc = acc
        preseq = seq
      end
      trace << "#{acc}\t#{preacc}\n"
      while 0 < trace.bytesize
        begin
          written = tracefp.write_nonblock(trace)
          trace = trace.byteslice(written..-1)
        rescue IO::WaitWritable
          if 1048576 < trace.bytesize
            IO.select(nil, [tracefp])
            retry
          end
        end
      end
    end
    infp.close
    tracefp.write(trace) if 0 < trace.bytesize
    tracefp.close
    outfp.write(out) if 0 < out.bytesize
    outfp.close
  end

  Process.waitpid(pid)
  sh "touch #{t.name}.trace"
end

file 'tmp/step2b.trace' => 'tmp/step2b'
