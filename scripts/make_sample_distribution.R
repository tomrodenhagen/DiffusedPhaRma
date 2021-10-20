x = rnorm(300,0,1)
k = ifelse(test = x <= 2, yes = "#8DD3C7", 
       no = ifelse(test = (x > 2 & x <= 58), yes =  "#FFFFB3", 
                 )          
hist(x,breaks=30,col=k)