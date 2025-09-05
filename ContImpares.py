def cont_impares(n):
    c,i = 0,1
    for i in range(n):
        if i % 2 == 1:
            c += 1
    print(c)

cont_impares(10)