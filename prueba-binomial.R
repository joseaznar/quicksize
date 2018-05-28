library("ggplot2")
library("grid")
library("gridExtra")

# primero ponemos los valores para correr las funciones
s <- 123456
set.seed(s)
iter <- 500
B <- 10
n0 <- 20

# creamos la funcion que genera una matriz de contingencia dado un tamagno de muestra
m.contingencia <- function(n,p1=0.1,p2=0.5){
  # para crear la matriz sacamos el numero de exitos t fracasos que se tienen con dos binomiales
  b1 <- rbinom(1,n,p1)
  b2 <- rbinom(1,n,p2)
  exitos <- c(b1,b2)
  rechazos <- c(n-b1,n-b2)
  datos <- c(exitos,rechazos)
  
  return(matrix(data=datos,nrow=2,ncol=2))
}

# el siguiente metodo simula y rechaza un muestra de dos binomiales independientes 
# dada una n usando la prueba exacta de fisher
sim.rechazo <- function(n){
  # primero obtenemos la matriz de contingencia
  cont <- m.contingencia(n)
  or <- 1.25
  # ahora calculamos la prueba exacta de fisher para esta matriz utilizando el ratio entre 0.1 y 0.5 e indicando que sabemos que
  # la p1 < p2
  fish <- fisher.test(x=cont, or=or)
  
  # .-
  if(fish$p.value<0.05){
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}

# el siguiente paso que haremos es hacer una funcion que nos diga el ratio de aceptacion que sera nuestra aproximacion 
# para el error de tipo dos
sim.beta <- function(B,n){
  cont <- 0
  for(i in 1:B){
    if(sim.rechazo(n)){
      cont <- cont + 1
    }
  }
  return(cont)
}

# finalmente creamos la funcion que obtiene el tamagno de la muestra dada una estimacion n0 inicial,un numero de iteraciones I,
# un numero de muestras por iteracion B, una beta objetivo y la delta para n a cada paso
quick.size <- function(n0,I,B,b,m){
  n <- n0
  N <- NULL
  n.aux <- NULL
  betas <- NULL
  llaves <- list()
  rechazos <- list()
  N <- list()
  # vamos a contar cuantos elementos se tienen
  num.elems <- 0
  
  # vamos a iterar el numero de veces propuesto
  for(i in 1:I){
    # metemos la n al arreglo
    n.aux[i] <- n
    
    # obtenemos  el numero de rechazos
    #rech <- test1(n,B)
    rech <- sim.beta(B,n)
    
    # revisamos si ya habiamos tenido esta n antes
    val <- match(n,llaves)
    if(is.na(val)){
      num.elems <- num.elems + 1
      
      # si no la habiamos tenido antes, la agregamos a las llaves y ponemos el valor correspondiente en las otras listas
      llaves[num.elems] <- n
      rechazos[num.elems] <- rech
      N[num.elems] <- B
      beta.aprox <- rech/B
    }
    else{
      # si ya lo teniamos, sumamos los valores correspondientes a el arreglo
      rechazos[val] <- as.numeric(rechazos[val][1]) + rech
      N[val] <- as.numeric(N[val]) + B
      
      # con estos nuevos valores calculamos la aproximacion para beta
      beta.aprox <- as.numeric(rechazos[val])/as.numeric(N[val])
    }
    betas[i] <- beta.aprox
    # ahora aumentamos, disminuimos o dejamos igual a n dependiendo del cociente obtenido con respecto al numero de rechazos
    if(beta.aprox < 1-b){
      n <- n + m
    }
    else{
      if(beta.aprox > 1-b){
        n <- n - m
      }
    }
  }
  return(list(n.aux,llaves,rechazos,N,betas))
}

# obtenemos el resultado deseado con las probabilidades preestablecidad
resp <- quick.size(n0,iter,B,0.2,1)

# generamos la grafica con los resultados para ver la convergencia
rew <- data.frame(matrix(c(unlist(resp[1]),1:iter),ncol=2))
colnames(rew) <- c("n","iteration")
grafica <- ggplot(data=rew, aes(iteration,n)) + geom_line() + geom_point(color="blue") + xlim(0, iter) + ggtitle(paste("B =",B," I = ", iter, " seed = ", s, " n0=", n0))

# ordenamos los resultados para mostrarlos como queremos
m <- matrix(c(unlist(resp[2]),unlist(resp[3]),unlist(resp[4]),as.numeric(unlist(resp[3]))/as.numeric(unlist(resp[4]))),ncol = 4)
m <- m[order(m[,1]),]
colnames(m) <- c("n","B*(n)", "B(n)", "P*(n)")

# imprimimos la grafica y la tabla juntos
grid.arrange(grafica,tableGrob(m))

