require 'levenshtein'

##
## Step 2a - exclusion of redundant reads by PCR & barcode matching
##

step2a_sources = [File.expand_path(CONF[LIBIDS[0]]['LAYOUT'])]
LIBIDS.each do |libid|
  CONF[libid]['FASTQS'].each_index do |runid|
    step2a_sources.unshift("tmp/#{libid}.#{runid}.step1b")
  end
end

file 'tmp/step2a' => step2a_sources do |t|
  pids = Array.new
  
  conf = CONF[LIBIDS[0]]
  umi = conf['UMI']
  barcode = conf['BARCODE']
  gap = conf['GAP']
  cdna = conf['CDNA']
  end5 = umi+barcode+gap+1
  end3 = umi+barcode+gap+cdna

  wells = Array.new
  barcodegaps = Array.new
  infp = open(t.sources[-1])
  while line = infp.gets
    well, barcodegap = line.rstrip.split(/\t/)
    wells.push(well)
    barcodegaps.push(barcodegap)
  end
  infp.close

  infifopaths = Array.new
  outfifopaths = Array.new
  PROCS.times do |i|
    infifopaths.push(mymkfifo('step2a-input-'))
    outfifopaths.push(mymkfifo('step2a-output-'))
  end

  pid = fork do
    fifos = Array.new
    infifopaths.each { |fifopath| fifos.push(open(fifopath, 'w')) }
    fifoidx = 0
    tracefp = open("| gsort --parallel=#{PROCS} -S #{PROCS}G | pigz -c > #{t.name}.trace", 'w')
    infp = open("| unpigz -c #{t.sources[0..-2].join(' ')} | gsort --parallel=#{PROCS} -S #{PROCS}G -k 4,4 -k 3,3r -m")
    line = infp.gets
    prelibid, tmpacc, preqv, preseq = line.rstrip.split(/\t/)
    preacc = "#{tmp=tmpacc.split(/:/); tmp[0..-2].join(':')}:#{end5}-#{end3}"
    tracefp.puts [preacc, preacc].join("\t")
    while line = infp.gets
      libid, tmpacc, qv, seq = line.rstrip.split(/\t/)
      acc = "#{tmp=tmpacc.split(/:/); tmp[0..-2].join(':')}:#{end5}-#{end3}"
      unless preseq == seq
        fifos[fifoidx].puts [prelibid, preacc, preqv, preseq].join("\t")
        fifoidx += 1
        fifoidx = 0 if fifoidx == PROCS
        prelibid = libid
        preacc = acc
        preqv = qv
        preseq = seq
      end
      tracefp.puts [preacc, acc].join("\t")
    end
    infp.close
    tracefp.close
    fifos[fifoidx % PROCS].puts [prelibid, preacc, preqv, preseq].join("\t")
    fifos.each { |fifo| fifo.close }
    Kernel.exit!
  end
  pids.push(pid)

  PROCS.times do |i|
    pid = fork do
      infp = open(infifopaths[i], 'r')
      outfp = open(outfifopaths[i], 'w')
      while line = infp.gets
        libid, acc, qv, seq = line.rstrip.split(/\t/)
        tmpdists = Hash.new
        barcodegaps.each_index do |idx|
          tmpdist = Levenshtein.distance(barcodegaps[idx], seq[umi, barcode+gap], threshold=2)
          dist = tmpdist.nil? ? 2 : tmpdist
          tmpdists[wells[idx]] = dist
          break if dist < 2
        end
        dists = tmpdists.sort { |a, b| a[1] <=> b[1] }
        well = dists[0][1] < 2 ? dists[0][0] : 'nobc'
        outfp.puts(["#{libid}.#{well}", acc, qv[end5-1, cdna], seq[end5-1, cdna]].join("\t"))
      end
      outfp.close
      infp.close
      Kernel.exit!
    end
    pids.push(pid)
  end

  outfp = open("| gsort --parallel=#{PROCS} -S #{PROCS}G -k 4,4 -k 3,3r | pigz -c > #{t.name}", 'w')
  fifos = Array.new
  outfifopaths.each { |fifopath| fifos.push(open(fifopath, 'r')) }
  fifodones = Hash.new
  fifos.cycle do |fifo|
    unless fifodones.key?(fifo)
      line = fifo.gets
      if line.nil?
        fifodones[fifo] = ''
        break if fifodones.keys.length == fifos.length
      else
        outfp.puts line.rstrip
      end
    end
  end
  fifos.each { |fifo| fifo.close }
  outfp.close

  pids.each do |pid|
    Process.waitpid(pid)
  end
end

task :clean_step2a do
  rm_rf "tmp/step2a"
  rm_rf "tmp/step2a.trace"
end
