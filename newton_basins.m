#funcao dada em outro arquivo (teremos que testar varias funcoes)
f = @myfunction;
#funcao derivada de "f". Dada em outro arquivo
df = @myderivative;
function newton basins (f, l, u, p)
    #{
    Acha as bacias de convergência da função f no dominio [l1, u1]×[l2, u2] e gera um arquivo output.txt que contém os 
    dados para a geração da imagem das bacias (pode usar gnuplot para gerar as imagens). 
    Os dados gerados preenchem uma imagem com p1 × p2 pixels. 
    #}
endfunction

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
    while()
    #invariante: x eh o ultimo numero calculado da sequencia do metodo de newton (i-esimo x da sequencia)
        x = x0 - (f(x0)/df(x0));
        i++;
        if (x - x0 <= epsilon || i >= iterations) 
        #se a diferenca entre o x e seu anterior for pequena o bastante ou ultrapassarmos o limite de iteraçoes do
        #algoritmo, devemos retornar x.
            return
        endif
        x0 = x;
    endwhile
    
    
endfunction