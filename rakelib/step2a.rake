##
## Step 2a - exclusion of redundant reads by PCR & barcode matching
##

require 'damerau-levenshtein'

def step2a_byLibraries_sources(path)
  libid = path.pathmap('%X').pathmap('%f')
  tmp = [getLayout(libid)]
  LIBRARIES[libid]['FASTQS'].each_index do |runid|
    tmp.push("tmp/#{libid}.#{runid}.step1b")
  end
  return tmp
end

rule '.step2a' => [->(path){ step2a_byLibraries_sources(path) }] do |t|
  umi = PREPROCESS['UMI']
  barcode = PREPROCESS['BARCODE']
  gap = PREPROCESS['GAP']
  cdna = PREPROCESS['CDNA']
  end5 = umi+barcode+gap+1
  end3 = umi+barcode+gap+cdna

  wells = Array.new
  barcodegaps = Array.new
  infp = open(t.source)
  while line = infp.gets
    well, barcodegap = line.rstrip.split(/\t/)
    wells.push(well)
    barcodegaps.push(barcodegap)
  end
  infp.close

  pid = fork do 
    outfp = open("| gsort --parallel=#{PROCS} -S 30% -t '\t' -k 4,4 -k 3,3r | pigz -c > #{t.name}", 'w')
    tracefp = open("| gsort --parallel=#{PROCS} -S 30% -t '\t' | pigz -c > #{t.name}.trace", 'w')
    infp = open("| unpigz -c #{t.sources[1..-1].join(' ')} | gsort --parallel=#{PROCS} -S 30% -t '\t' -k 4,4 -k 3,3r -m")
    preseq = ''
    prebcseq = ''
    preacc = ''
    while line = infp.gets
      ## puts "#{infp.lineno} #{`date`}" if infp.lineno % 100000 == 0
      libid, tmpacc, qv, seq = line.rstrip.split(/\t/)
      acc = "#{tmp=tmpacc.split(/:/); tmp[0..-2].join(':')}:#{end5}-#{end3}"
      if preseq != seq
        bcseq = seq[umi, barcode+gap]
        if prebcseq != bcseq
          bcidx = nil
          barcodegaps.each_index do |idx|
            dist = DamerauLevenshtein.distance(barcodegaps[idx], bcseq, 0, max_distance=2)
            if !dist.nil? && dist < 2
              bcidx = idx
              break
            end
          end
          prebcseq = bcseq
        end
        begin
          outfp.write_nonblock("#{libid}.#{bcidx.nil? ? 'nobc' : wells[bcidx]}\t#{acc}\t#{qv[end5-1, cdna]}\t#{seq[end5-1, cdna]}\n")
        rescue IO::WaitWritable
          IO.select(nil, [outfp])
          retry
        end
        preacc = acc
        preseq = seq
      end
      begin
        tracefp.write_nonblock("#{preacc}\t#{acc}\n")
      rescue IO::WaitWritable
        IO.select(nil, [tracefp])
        retry
      end
    end
    infp.close
    tracefp.close
    outfp.close
  end

  Process.waitpid(pid)
  sh "touch #{t.name}.trace"
end

rule /\.step2a\.trace$/ => [->(path){ path.pathmap("%X") }]
