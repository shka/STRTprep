##
## Step 3f - Reads and the normalized reads in the qualified samples
##

def step3f_job(t, colname)
  qualified = Hash.new
  infp = open("| grep -v ^LIBRARY #{t.source}", 'rt')
  while line = infp.gets
    cols = line.rstrip.split(/,/)
    qualified["#{cols[0]}.#{cols[1]}"] = ''
  end
  infp.close

  header = [true]
  outfp = open("| gzip -c > #{t.name}", 'w')
  infp = open("| gunzip -c #{t.sources[1]}")
  while line = infp.gets
    cols = line.rstrip.split(/\t/)
    if header.length == 1
      tmp = [colname]
      cols[1..-1].each do |colname|
        libwellid, name = colname.split(/\|/)
        header.push(qualified.key?(libwellid))
        tmp.push(colname) if qualified.key?(libwellid)
      end
      outfp.puts tmp.join("\t")
    else
      tmp = [cols[0]]
      sum = 0
      nonzero = 0
      1.upto(header.length-1) do |i|
        if header[i]
          cnt = cols[i].to_i
          tmp.push(cnt)
          sum += cnt
          nonzero += 1 if cnt > 0
        end
      end
      if colname == 'Gene'
        outfp.puts tmp.join("\t") if sum > 0
      else
        outfp.puts tmp.join("\t") if sum > 0 && nonzero > 2
      end
    end
  end
  infp.close
  outfp.close
end

file 'out/byGene/reads.txt.gz' =>
     ['out/byGene/samples.csv', 'out/byGene/reads_all.txt.gz'] do |t|
  step3f_job(t, 'Gene')
end

rule /\/reads\.RData$/ => [->(path){ path.sub(/\.RData$/, '.txt.gz') }] do |t|
  sh <<EOF
echo "tmp.reads <- as.matrix(read.table('#{t.source}', header=T, sep='\\t', quote='', row.names=1, check.names=F)); reads <- tmp.reads[order(rownames(tmp.reads)), order(colnames(tmp.reads))]; save(reads, file='#{t.name}', compress='gzip')"\
| R --vanilla --quiet > #{t.name}.log 2>&1
EOF
end

rule /\/nreads\.txt\.gz$/ => [->(path){ path.sub(/\/nreads\.txt\.gz$/, '/reads.RData') }] do |t|
  sh "R --vanilla --quiet --args #{t.source} #{t.name} < bin/_step3f_normalization.R > #{t.name}.log 2>&1"
end

rule /\/nreads\.RData$/ => [->(path){ path.sub(/\.RData$/, '.txt.gz') }] do |t|
  sh <<EOF
echo "nreads <- as.matrix(read.table('#{t.source}', header=T, sep='\\t', quote='', row.names=1, check.names=F)); save(nreads, file='#{t.name}', compress='gzip')"\
| R --vanilla --quiet > #{t.name}.log 2>&1
EOF
end
