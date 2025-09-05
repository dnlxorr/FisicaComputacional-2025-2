def es_primo(n):
  if n < 2:
    return False
  for i in range(2, int(n**0.5) + 1):
    if n % i == 0:
          return False
  return True

contador = 0
for x in range(1, 501):  #Desde 1 hasta 500
  if es_primo(x):
    contador += 1

print("Cantidad de números primos entre 1 y 500:", contador)