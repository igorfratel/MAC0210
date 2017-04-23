1;
function newton_basins (f, df, l, u, p)
    #{
    Acha as bacias de convergência da função f no dominio [l1, u1]×[l2, u2] e gera um arquivo output.txt que contém os 
    dados para a geração da imagem das bacias (pode usar gnuplot para gerar as imagens). 
    Os dados gerados preenchem uma imagem com p1 × p2 pixels. 
    #}
    graphics_toolkit ("gnuplot");
    axis_config = [l{1,1}, u{1,1}, l{1,2}, u{1,2}] ;
    axis(axis_config);
    epsilon = 0.01;
    iterations = 100;
    vertical = (u{1,2} - l{1,2})/p{1,2};
    horizontal = (u{1,1} - l{1,1})/p{1,1};
    color_matrix = cell(p{1,1}, p{1,2});
    y = l{1,2};
    j = p{1,1};
    while (y <= u{1,2})
        x = l{1,1} + y;
        i = 1;
        while (x <= u{1,1} + y)
            result = newton(f, df, x, epsilon, iterations);
            #ver se o valor de result eh alguma raiz (tem q ter um vetor de raizes)
            color_matrix{j}{i} = #????
            x += horizontal;
            i++;
        endwhile
        j--;
        y += vertical;
    endwhile
    #normalizar a matriz;
    #gerar arquivo texto;
    imagesc(color_matrix);
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
#funcao dada em outro arquivo (teremos que testar varias funcoes)
f = @myfunction;
df = @myderivative;
l = {-2, -2*i};
u = {2, 2*i};
p = {2,2};
newton_basins(f, df, l, u, p);