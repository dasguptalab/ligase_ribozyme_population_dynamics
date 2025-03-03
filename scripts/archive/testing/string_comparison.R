numMatch <- mapply(
  function(x, y) {
    len <- length(x)
    sum(x[1:len] == y[1:len])
  }, 
  strsplit("GAATGCTGCCAACCGTGCGGGCTAATTGGCAGACTGAGCT", ''), 
  strsplit("TGATGAATCGGCATACGTGGGTCAGAGTCATAGTGCGACA", '')
)

40-numMatch

