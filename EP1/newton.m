function x = newton (f, df, x0, epsilon, iterations)
    #{
    Aplica o método de Newton para achar uma raiz da função f (com primeira derivada df), partindo do ponto x0.
    Os parâmetros f e df devem ser apontadores a funções em Octave que implementem f e sua derivada primeira, 
    respectivamente. 
    A função que implementa o método de Newton pode ter parâmetros adicionais, relacionados com tolerâncias para o 
    critério de parada ou um máximo de iterações para o caso em que a convergência não ocorra.
    #}
    
    #Adicionar criterio para quando df(x) == 0!!!
    i = 0;
    while(1)
    #invariante: x eh o ultimo numero calculado da sequencia do metodo de newton (i-esimo x da sequencia)
        if (df(x0) != 0)
            x = x0 - (f(x0)/df(x0));
            i++;
        else
            x0
            return 
        endif
        if (x - x0 <= epsilon || i >= iterations) 
        #se a diferenca entre o x e seu anterior for pequena o bastante ou ultrapassarmos o limite de iteraçoes do
        #algoritmo, devemos retornar x.
            x
            return
        endif
        x0 = x;
    endwhile
    x  
endfunction