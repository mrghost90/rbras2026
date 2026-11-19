# Função para a construção de Látice quadrado
# em que k é número primo:
Latice_Primo <- function(k){
  t <- k^2; r <- k+1; b <- r*k;  l <- 1
  # Conjunto básico A:
  A_X <- t(matrix(1:t,k,k)) 
  A_Y <- t(A_X) 
  A_Z <- A_Y
  A_Z_base <- A_Y
  planos <- rbind(A_X,A_Y)
  jj <- NULL
  for(count in 3:r){
    A_Z <- A_Z_base
    for(i in 2:k){#colunas
      for(j in 1:k){#linhas
        jj <- (j+i-1)%%k
        if(jj==0){
          jj<-k
        }
        A_Z[j,i] <- A_Z_base[jj,i]
      }
    }
    A_Z_base <- A_Z
    planos <- rbind(planos,A_Z)
  }
  design <- t(planos)
  B  <- factor(kronecker(1:b,rep(1,k)))
  Tr <- factor(c(design))
  
  R <- kronecker(1:(k+1),rep(1,k^2))
  planilha <- data.frame(R,B,Tr)
  Z <- model.matrix(~ -1+Tr)
  X <- model.matrix(~ -1+factor(B))
  
  XZ <- cbind(X,Z)
  ZXXZ <- t(XZ)%*%XZ
  
  N <- t(X)%*%Z
  LAMBDA <- t(N)%*%N
  if(sum(LAMBDA==l)!=(t^2-t)){
    print("CUIDADO - DELINEAMENTO NÃO É BIB!")
    print("Tem certeza que k é primo?")
  }else{
    return(list(planilha=planilha,ZXXZ=ZXXZ,LAMBDA=LAMBDA))
  }
}

k <- 11
saida <- Latice_Primo(k)

saida$LAMBDA

write(colnames(saida$planilha),file="planilha_Miguel.csv",ncol=3)
write(t(saida$planilha),file="planilha_Miguel.csv",ncol=3,append=TRUE)
