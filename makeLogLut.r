
limit <- 0.003

tableSize <- 10000000
p <- seq(0,limit,by=1/tableSize)

cat(sprintf("#define ENTROPY_TABLE_SIZE %f\n",tableSize))

cat(sprintf("#define ENTROPY(_p) ( (_p) > %f ?  -(_p)*log2(_p) : entropyLUT[(int)((double)(_p)*ENTROPY_TABLE_SIZE)])\n",limit))

cat(sprintf("double entropyLUT[] = { %f", 0 ))

for(i in 2:length(p))
    cat(sprintf(",%10.7f", - p[i] * log2(p[i])))
cat("};")
