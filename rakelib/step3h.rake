##
## Step 3h - differential expression analysis
##

##

def step3h_testN_sources(path)
  return ["out/byGene/samples.csv",
          "#{path.pathmap('%d')}/reads.RData",
          "#{path.pathmap('%d')}/nreads.RData"]
end

rule /\/sources\d+\.RData$/ => [->(p){ step3h_testN_sources(p) }] do |t|
  idx = /sources(\d+)\.RData$/.match(t.name).to_a[1]
  sh "bin/_step3h_extract_sources.R #{idx} #{t.sources.join(' ')} #{t.name} > #{t.name}.log 2>&1"
end

##

def step3h_diffexpN_sources(path)
  root = path.pathmap('%d')
  idx = /diffexp(\d+)\.txt\.gz$/.match(path).to_a[1]
  return ["#{root}/sources#{idx}.RData"]
end

require 'csv'

rule /\/diffexp\d+\.txt\.gz$/ =>
                         [->(path){ step3h_diffexpN_sources(path) }] do |t|
  idx = /diffexp(\d+)\.txt\.gz$/.match(t.name).to_a[1]
  tmp = t.source.sub(/^out/, 'tmp')
  if (!File.exist?(t.name) ||
      !File.exist?(tmp) ||
      `gmd5sum #{t.source} #{tmp} | gcut -d ' ' -f 1 | guniq | gwc -l`.to_i != 1)
    sh "bin/_step3h_test_diffexp.R #{idx} #{t.source} #{t.name} #{t.name.pathmap('%d')} > #{t.name}.log 2>&1"
    sh "gcp --preserve=timestamps #{t.source} #{tmp}"
  else
    puts "... skipped #{t.name}, since the target samples were identical with the previous run"
    sh "touch #{t.name}"
  end
end

##

rule /\/fluctuation\d+\.txt\.gz$/ =>
                             [->(p){ p.sub('/fluctuation', '/diffexp') }]

##

step3h_sources = ['out/byGene/reads.txt.gz',
                  'out/byGene/nreads.txt.gz',
                  'out/byGene/fluctuation.txt.gz']
begin
  infp = open('src/samples.csv', 'rt')
  colnames = infp.gets.rstrip.split(',')
  tmp = colnames.select { |colname| /^CLASS\.\d+$/.match(colname) }
  classes = tmp.map { |cls| /^CLASS\.(\d+)$/.match(cls).to_a[1] }
  infp.close
  classes.each do |cls|
    step3h_sources.push("out/byGene/diffexp#{cls}.txt.gz")
  end
end

def step3h_job(t)
  ofs = /regions/ =~ t.sources[3] ? 5 : 3
  
  annotation = nil
  if ofs == 5
    annotation = Hash.new
    
    infp = open("| unpigz -c #{t.sources[3]}")
    while line = infp.gets
      cols = line.rstrip.split(/\t/)
      start = cols[1].to_i+1
      annotation[cols[3]] = ["#{cols[0]}:#{start}-#{cols[2]}", cols[5], cols[0], start, cols[2]]
    end
    infp.close
    
    infp = open("| unpigz -c #{t.sources[4]}")
    while line = infp.gets
      cols = line.rstrip.split(/\t/)
      annotation[cols[0]] = [cols[4], cols[6]] + annotation[cols[0]] + [cols[2], cols[5]]
    end
    infp.close
  end
  
  diffexps = Hash.new
  header_diffexps = Array.new
  t.sources[ofs..-1].each_index do |idx|
    tmp = Hash.new
    infp = open("| unpigz -c #{t.sources[idx+ofs]}")
    clsnum = /\/diffexp(\d+)/.match(t.sources[idx+ofs]).to_a[1].to_i
    header = infp.gets
    while line = infp.gets
      cols = line.rstrip.split(/\t/)
      tmp[cols[0]] = cols[1..-1].join(',')
    end
    infp.close
    infp = open("| unpigz -c #{t.sources[idx+ofs].sub('/diffexp', '/fluctuation')}")
    header = infp.gets
    while line = infp.gets
      cols = line.rstrip.split(/\t/)
      tmp[cols[0]] = [tmp[cols[0]]] + (cols[1] != 'NA' ? [cols[1], cols[2]] : [1, 'NA']) if tmp.key?(cols[0])
    end
    infp.close
    diffexps[idx] = tmp
    header_diffexps.push("diffexpScore.#{clsnum}")
    header_diffexps.push("pvalue.#{clsnum}")
    header_diffexps.push("qvalue.#{clsnum}")
    header_diffexps.push("fluctuation.#{clsnum}")
    header_diffexps.push("fluctuationScore.#{clsnum}")
  end

  fluctuation = Hash.new
  infp = open("| unpigz -c #{t.sources[2]}")
  header_fluctuation = infp.gets
  while line = infp.gets
    gene, pvalue, score = line.rstrip.split(/\t/)
    fluctuation[gene] = pvalue != 'NA' ? [pvalue, score] : [1, 'NA']
  end
  infp.close
  
  nreads = Hash.new
  infp = open("| unpigz -c #{t.sources[1]}")
  header_nreads = infp.gets.rstrip.split(/\t/)
  1.upto(header_nreads.length-1) do |i|
    header_nreads[i] = "N|#{header_nreads[i]}"
  end
  while line = infp.gets
    cols = line.rstrip.split(/\t/)
    nreads[cols[0]] = cols[1..-1]
  end
  infp.close

  outfp = open(t.name, 'w')
  infp = open("| unpigz -c #{t.source}")
  header_reads = infp.gets.rstrip.split(/\t/)
  1.upto(header_reads.length-1) do |i|
    header_reads[i] = "R|#{header_reads[i]}"
  end
  outfp.puts ([header_reads[0]] + (annotation.nil? ? [] : ['Gene', 'Location', 'Region', 'Str', 'Chr', 'Start', 'Stop', 'Peak', 'Transcript']) + ['fluctuation.global', 'fluctuationScore.global'] + header_diffexps + header_nreads[1..-1] + header_reads[1..-1]).join(',')
  while line = infp.gets
    cols = line.rstrip.split(/\t/)
    row_diffexps = Array.new
    diffexps.keys.sort.each do |idx|
      tmp = diffexps[idx]
      row_diffexps.push(tmp.key?(cols[0]) ? tmp[cols[0]] : ",,,,")
    end
    if annotation.nil?
      outfp.puts ([cols[0]] + fluctuation[cols[0]] + row_diffexps + nreads[cols[0]] + cols[1..-1]).join(',')
    elsif annotation[cols[0]].length == 9
      outfp.puts ([cols[0]] + annotation[cols[0]] + fluctuation[cols[0]] + row_diffexps + nreads[cols[0]] + cols[1..-1]).join(',')
    end
  end
  infp.close
  outfp.close
end

file 'out/byGene/diffexp.csv' => step3h_sources do |t|
  step3h_job(t)
end
