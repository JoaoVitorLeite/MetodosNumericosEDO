


// LaTeX $\frac{dy}{dx} = \frac{3}{4} - \frac{3 \cdot y}{100}$
//Função do PVI
function r = f_pvi(x,y)
    r = 3/4 - (3*y)/100
    //r = y-x
    //r = x - y + 2
endfunction



// LaTeX $\frac{dy}{dx} = 25 + 25 \cdot e^{-3 \cdot x / 100}$
//Função resolução da EDO
function r = f_solucao(x)
    r = 25 + 25*(%e^((-3*x)/100))
    //r = %e^x + x + 1
    //r = %e^(-x) + x + 1
endfunction



// LaTeX $x_{n+1} = x_{n} + h$
// LaTeX $y_{n+1} = y_{n} + h \cdot f(x_n, y_n)$
function [x, y, tempo] = metodo_euler(a, b, N, h, x0, y0, log_output)
    
    disp("METODO DE EULER")
    
    if ~exists('log_output', 'local') then
        log_output = %T
    end
    
    if log_output then
        [fd, err] = mopen("log_euler.txt", "a+")
        if err ~= 0 then
            return
        end
    end
    
    if log_output then
        mfprintf(fd, "%f,%f,%f,%f,%f,%f\n", a,b,N,h,x0,y0)
        mfprintf(fd, "x_i,y_i_aproximado,y_i_analitico,erro\n")
    end
    
    tic()
    
    x = zeros(N+1)
    y = zeros(N+1)
    y_analitico = zeros(N+1)
    
    x(1) = x0
    y(1) = y0
    y_analitico(1) = y0
        
    for i = 1:N
        x(i+1) = x(i) + h
        y(i+1) = y(i) + h*f_pvi(x(i), y(i))
        y_analitico(i+1) = f_solucao(x(i+1))
        erro = abs(y_analitico(i+1) - y(i+1))
        
        disp([x(i+1) y(i+1) y_analitico(i+1) erro])
        
        if log_output then
            mfprintf(fd, "%.16f,%.16f,%.16f,%.16f\n", x(i+1), y(i+1), y_analitico(i+1), erro)
        end
        
    end
    
    tempo = toc()
    
    if log_output then
        mfprintf(fd, "Tempo=%f\n\n", tempo)
        mclose(fd);
    end
    
    clf
    plot2d(x', [y y_analitico], style=[20 1], leg="PVI@ANALITICO")
    title("SOLUCAO", "color", "blue")
    xlabel("x", "color", "blue")
    ylabel("y", "color", "blue")
    
endfunction



// LaTeX $x_{n+1} = x_{n} + h$
// LaTeX $K_1 = f(x_n, y_n)$
// LaTeX $y_{n+1} = y_n + h \cdot K_1$
function [x, y, tempo] = metodo_runge_kutta_1_ordem(a, b, N, h, x0, y0, log_output)
    
    disp("METODO DE RUNGE KUTTA DE ORDEM 1")
    
    if ~exists('log_output', 'local') then
        log_output = %T
    end
    
    if log_output then
        [fd, err] = mopen("log_rk1.txt", "a+")
        if err ~= 0 then
            return
        end
    end
    
    if log_output then
        mfprintf(fd, "%f,%f,%f,%f,%f,%f\n", a,b,N,h,x0,y0)
        mfprintf(fd, "x_i,y_i_aproximado,y_i_analitico,erro\n")
    end
    
    tic()
    
    x = zeros(N+1)
    y = zeros(N+1)
    y_analitico = zeros(N+1)
    
    x(1) = x0
    y(1) = y0
    y_analitico(1) = y0
    
    for i = 1:N
        
        x(i+1) = x(i) + h
        k1 = f_pvi(x(i), y(i))
        y(i+1) = y(i) + h*k1
        y_analitico(i+1) = f_solucao(x(i+1))
        erro = abs(y_analitico(i+1) - y(i+1))
        
        disp([x(i+1) y(i+1) y_analitico(i+1) erro])
        
        if log_output then
            mfprintf(fd, "%.16f,%.16f,%.16f,%.16f\n", x(i+1), y(i+1), y_analitico(i+1), erro)
        end
    end
    
    tempo = toc()
    
    if log_output then
        mfprintf(fd, "Tempo=%f\n\n", tempo)
        mclose(fd);
    end
    
    clf
    plot2d(x', [y y_analitico], style=[20 1], leg="PVI@ANALITICO")
    title("SOLUCAO", "color", "blue")
    xlabel("x", "color", "blue")
    ylabel("y", "color", "blue") 
    
endfunction



// LaTeX $x_{n+1} = x_{n} + h$
// LaTeX $K_1 = f(x_n, y_n)$
// LaTeX $K_2 = f(x_n + h, y_n + h \cdot K_1)$
// LaTeX $y_{n+1} = y_n + \frac{h}{2} \cdot (K_1 + K_2)$
function [x, y, tempo] = metodo_runge_kutta_2_ordem(a, b, N, h, x0, y0, log_output)
    
    disp("METODO DE RUNGE KUTTA DE ORDEM 2")
    
    if ~exists('log_output', 'local') then
        log_output = %T
    end
    
    if log_output then
        [fd, err] = mopen("log_rk2.txt", "a+")
        if err ~= 0 then
            return
        end
    end
    
    if log_output then
        mfprintf(fd, "%f,%f,%f,%f,%f,%f\n", a,b,N,h,x0,y0)
        mfprintf(fd, "x_i,y_i_aproximado,y_i_analitico,erro\n")
    end
    
    tic()
    
    x = zeros(N+1)
    y = zeros(N+1)
    y_analitico = zeros(N+1)
    
    x(1) = x0
    y(1) = y0
    y_analitico(1) = y0
    
    for i = 1:N
        
        x(i+1) = x(i) + h
        k1 = f_pvi(x(i), y(i))
        k2 = f_pvi(x(i) + h, y(i) + h*k1)
        y(i+1) = y(i) + (h/2)*(k1 + k2)
        y_analitico(i+1) = f_solucao(x(i+1))
        erro = abs(y_analitico(i+1) - y(i+1))
        
        disp([x(i+1) y(i+1) y_analitico(i+1) erro])
        
        if log_output then
            mfprintf(fd, "%.16f,%.16f,%.16f,%.16f\n", x(i+1), y(i+1), y_analitico(i+1), erro)
        end
        
    end
    
    tempo = toc()
    
    if log_output then
        mfprintf(fd, "Tempo=%f\n\n", tempo)
        mclose(fd);
    end
    
    clf
    plot2d(x', [y y_analitico], style=[20 1], leg="PVI@ANALITICO")
    title("SOLUCAO", "color", "blue")
    xlabel("x", "color", "blue")
    ylabel("y", "color", "blue")
    
endfunction



// LaTeX $x_{n+1} = x_{n} + h$
// LaTeX $K_1 = f(x_n, y_n)$
// LaTeX $K_2 = f(x_n + \frac{h}{2}, y_n + \frac{h}{2} \cdot K_1)$
// LaTeX $K_3 = f(x_n + h, y_n + 2 \cdot h \cdot K_2 - h \cdot K_1)$
// LaTeX $y_{n+1} = y_n + \frac{h}{6} \cdot (K_1 + 4 \cdot K_2 + K_3)$
function [x, y, tempo] = metodo_runge_kutta_3_ordem(a, b, N, h, x0, y0, log_output)
    
    disp("METODO DE RUNGE KUTTA DE ORDEM 3")
    
    if ~exists('log_output', 'local') then
        log_output = %T
    end
    
    if log_output then
        [fd, err] = mopen("log_rk3.txt", "a+")
        if err ~= 0 then
            return
        end
    end
    
    if log_output then
        mfprintf(fd, "%f,%f,%f,%f,%f,%f\n", a,b,N,h,x0,y0)
        mfprintf(fd, "x_i,y_i_aproximado,y_i_analitico,erro\n")
    end
    
    tic()
    
    x = zeros(N+1)
    y = zeros(N+1)
    y_analitico = zeros(N+1)
    
    x(1) = x0
    y(1) = y0
    y_analitico(1) = y0
    
    for i = 1:N
        
        x(i+1) = x(i) + h
        k1 = f_pvi(x(i), y(i))
        k2 = f_pvi(x(i) + (h/2), y(i) + (h/2)*k1)
        k3 = f_pvi(x(i) + h, y(i) + 2*h*k2 - h*k1)
        y(i+1) = y(i) + (h/6)*(k1 + 4*k2 + k3)
        y_analitico(i+1) = f_solucao(x(i+1))
        erro = abs(y_analitico(i+1) - y(i+1))
        
        disp([x(i+1) y(i+1) y_analitico(i+1) erro])
        
        if log_output then
            mfprintf(fd, "%.16f,%.16f,%.16f,%.16f\n", x(i+1), y(i+1), y_analitico(i+1), erro)
        end
        
    end
    
    tempo = toc()
    
    if log_output then
        mfprintf(fd, "Tempo=%f\n", tempo)
        mclose(fd);
    end
    
    clf
    plot2d(x', [y y_analitico], style=[20 1], leg="PVI@ANALITICO")
    title("SOLUCAO", "color", "blue")
    xlabel("x", "color", "blue")
    ylabel("y", "color", "blue")
    
endfunction



// LaTeX $x_{n+1} = x_{n} + h$
// LaTeX $K_1 = f(x_n, y_n)$
// LaTeX $K_2 = f(x_n + \frac{h}{2}, y_n + \frac{h}{2} \cdot K_1)$
// LaTeX $K_3 = f(x_n + \frac{h}{2}, y_n + \frac{h}{2} \cdot K_2)$
// LaTeX $y_{n+1} = y_n + \frac{h}{6} \cdot (K_1 + 2 \cdot K_2 + 2 \cdot K_3 + K_4)$
function [x, y, tempo] = metodo_runge_kutta_4_ordem(a, b, N, h, x0, y0, log_output)
    
    disp("METODO DE RUNGE KUTTA DE ORDEM 4")
    
    if ~exists('log_output', 'local') then
        log_output = %T
    end
    
    if log_output then
        [fd, err] = mopen("log_rk4.txt", "a+")
        if err ~= 0 then
            return
        end
    end
    
    if log_output then
        mfprintf(fd, "%f,%f,%f,%f,%f,%f\n", a,b,N,h,x0,y0)
        mfprintf(fd, "x_i,y_i_aproximado,y_i_analitico,erro\n")
    end
    
    tic()
    
    x = zeros(N+1)
    y = zeros(N+1)
    y_analitico = zeros(N+1)
    
    x(1) = x0
    y(1) = y0
    y_analitico(1) = y0
    
    for i = 1:N
        
        x(i+1) = x(i) + h
        k1 = f_pvi(x(i), y(i))
        k2 = f_pvi(x(i) + (h/2), y(i) + (h/2)*k1)
        k3 = f_pvi(x(i) + (h/2), y(i) + (h/2)*k2)
        k4 = f_pvi(x(i) + h, y(i) + h*k3)
        y(i+1) = y(i) + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
        y_analitico(i+1) = f_solucao(x(i+1))
        erro = abs(y_analitico(i+1) - y(i+1))
        
        disp([x(i+1) y(i+1) y_analitico(i+1) erro])
        
        if log_output then
            mfprintf(fd, "%.16f,%.16f,%.16f,%.16f\n", x(i+1), y(i+1), y_analitico(i+1), erro)
        end
        
    end
    
    tempo = toc()
    
    if log_output then
        mfprintf(fd, "Tempo=%f\n\n", tempo)
        mclose(fd);
    end
    
    clf
    plot2d(x', [y y_analitico], style=[20 1], leg="PVI@ANALITICO")
    title("SOLUCAO", "color", "blue")
    xlabel("x", "color", "blue")
    ylabel("y", "color", "blue")
    
endfunction



function [x, y, tempo] = metodo_runge_kutta(p, a, b, N, h, x0, y0, log_output)
    
    if ~exists('log_output', 'local') then
        log_output = %T
    end
    
    select p
    case 1 then
        [x, y, tempo] = metodo_runge_kutta_1_ordem(a, b, N, h, x0, y0, log_output)
    case 2 then
        [x, y, tempo] = metodo_runge_kutta_2_ordem(a, b, N, h, x0, y0, log_output)
    case 3 then
        [x, y, tempo] = metodo_runge_kutta_3_ordem(a, b, N, h, x0, y0, log_output)
    case 4 then
        [x, y, tempo] = metodo_runge_kutta_4_ordem(a, b, N, h, x0, y0, log_output)
    else
        break        
    end
    
endfunction

//TESTANDO TUDO ANTERIORMENTE


function [x, y, tempo] = metodo_passo_unico(a, b, N, h, x0, y0, method)
    
    select method
    case "EULER" then
        [x, y, tempo] = metodo_euler(a, b, N, h, x0, y0, %F)
    case "RK1" then
        [x, y, tempo] = metodo_runge_kutta_1_ordem(a, b, N, h, x0, y0, %F)
    case "RK2" then
        [x, y, tempo] = metodo_runge_kutta_2_ordem(a, b, N, h, x0, y0, %F)
    case "RK3" then
        [x, y, tempo] = metodo_runge_kutta_3_ordem(a, b, N, h, x0, y0, %F)
    case "RK4" then
        [x, y, tempo] = metodo_runge_kutta_4_ordem(a, b, N, h, x0, y0, %F)
    else
        break
    end
    
endfunction



// LaTeX $x_{n+1} = x_{n} + h$
// LaTeX $y_{n+1} = y_n + \frac{h}{2} \cdot (3 \cdot f(x_n, y_n) - f(x_{n-1}, y_{n-1})$
function [x,y,tempo] = metodo_adams_bashforth_passo_2(a, b, N, h, x0, y0, method, log_output)
    
    disp("METODO DE ADAMS BASHFORTH DE PASSO 2")
    
    if ~exists('log_output', 'local') then
        log_output = %T
    end
    
    if log_output then
        [fd, err] = mopen("log_adams_bashforth_passo_2.txt", "a+")
        if err ~= 0 then
            return
        end
    end
    
    if log_output then
        mfprintf(fd, "%f,%f,%f,%f,%f,%f,%s\n", a,b,N,h,x0,y0,method)
        mfprintf(fd, "x_i,y_i_aproximado,y_i_analitico,erro\n")
    end
    
    tic()
    
    disp("METODO DE PASSO UNICO")
    
    [x_aux, y_aux, tempo_aux] = metodo_passo_unico(a, b, 1, h, x0, y0, method)
        
    disp("METODO DE PASSO MULTIPLO")
    
    x = zeros(N+1)
    y = zeros(N+1)
    y_analitico = zeros(N+1)
    x(1:2) = x_aux(1:2)
    y(1:2) = y_aux(1:2)
    y_analitico(1:2) = f_solucao(x_aux(1:2))
    
    if log_output then
        for i = 2:2
            erro = abs(y_analitico(i) - y(i))
            mfprintf(fd, "%.16f,%.16f,%.16f,%.16f\n", x(i), y(i), y_analitico(i), erro)
        end
    end 
    
    for i = 2:N
        
        x(i+1) = x(i) + h
        y(i+1) = y(i) + (h/2)*(3*f_pvi(x(i), y(i)) - f_pvi(x(i-1), y(i-1)))
        y_analitico(i+1) = f_solucao(x(i+1))
        erro = abs(y_analitico(i+1) - y(i+1))
        
        disp([x(i+1) y(i+1) y_analitico(i+1) erro])
        if log_output then
            mfprintf(fd, "%.16f,%.16f,%.16f,%.16f\n", x(i+1), y(i+1), y_analitico(i+1), erro)
        end
        
    end
    
    tempo = toc()
    
    if log_output then
        mfprintf(fd, "Tempo=%f\n\n", tempo)
        mclose(fd);
    end
    
    
    clf
    plot2d(x', [y y_analitico], style=[20 1], leg="PVI@ANALITICO")
    title("SOLUCAO", "color", "blue")
    xlabel("x", "color", "blue")
    ylabel("y", "color", "blue")
    
endfunction



// LaTeX $x_{n+1} = x_{n} + h$
// LaTeX $y_{n+1} = y_n + \frac{h}{12} \cdot (23 \cdot f(x_n, y_n) - 16 \cdot f(x_{n-1}, y_{n-1}) + 5 \cdot f(x_{n-2}, y_{n-2}))$
function [x,y,tempo] = metodo_adams_bashforth_passo_3(a, b, N, h, x0, y0, method, log_output)
    
    disp("METODO DE ADAMS BASHFORTH DE PASSO 3")
    
    if ~exists('log_output', 'local') then
        log_output = %T
    end
    
    if log_output then
        [fd, err] = mopen("log_adams_bashforth_passo_3.txt", "a+")
        if err ~= 0 then
            return
        end
    end
    
    if log_output then
        mfprintf(fd, "%f,%f,%f,%f,%f,%f,%s\n", a,b,N,h,x0,y0,method)
        mfprintf(fd, "x_i,y_i_aproximado,y_i_analitico,erro\n")
    end
    
    tic()
    
    disp("METODO DE PASSO UNICO")
    
    [x_aux, y_aux, tempo_aux] = metodo_passo_unico(a, b, 2, h, x0, y0, method)
    
    disp("METODO DE PASSO MULTIPLO")
    
    x = zeros(N+1)
    y = zeros(N+1)
    y_analitico = zeros(N+1)
    x(1:3) = x_aux(1:3)
    y(1:3) = y_aux(1:3)
    y_analitico(1:3) = f_solucao(x_aux(1:3)) 
    
    if log_output then
        for i = 2:3
            erro = abs(y_analitico(i) - y(i))
            mfprintf(fd, "%.16f,%.16f,%.16f,%.16f\n", x(i), y(i), y_analitico(i), erro)
        end
    end   
    
    for i = 3:N
        
        x(i+1) = x(i) + h
        y(i+1) = y(i) + (h/12)*(23*f_pvi(x(i), y(i)) - 16*f_pvi(x(i-1), y(i-1)) + 5*f_pvi(x(i-2), y(i-2)))
        y_analitico(i+1) = f_solucao(x(i+1))
        erro = abs(y_analitico(i+1) - y(i+1))
        
        disp([x(i+1) y(i+1) y_analitico(i+1) erro])
        if log_output then
            mfprintf(fd, "%.16f,%.16f,%.16f,%.16f\n", x(i+1), y(i+1), y_analitico(i+1), erro)
        end
        
    end
    
    tempo = toc()
    
    if log_output then
        mfprintf(fd, "Tempo=%f\n\n", tempo)
        mclose(fd);
    end
    
    clf
    plot2d(x', [y y_analitico], style=[20 1], leg="PVI@ANALITICO")
    title("SOLUCAO", "color", "blue")
    xlabel("x", "color", "blue")
    ylabel("y", "color", "blue")
    
endfunction



// LaTeX $x_{n+1} = x_{n} + h$
// LaTeX $y_{n+1} = y_n + \frac{h}{24} \cdot (55 \cdot f(x_n, y_n) - 59 \cdot f(x_{n-1}, y_{n-1}) + 37 \cdot f(x_{n-2}, y_{n-2})            - 9 \cdot f(x_{n-3}, y_{n-3}))$
function [x,y,tempo] = metodo_adams_bashforth_passo_4(a, b, N, h, x0, y0, method, log_output)
    
    disp("METODO DE ADAMS BASHFORTH DE PASSO 4")
    
    if ~exists('log_output', 'local') then
        log_output = %T
    end
    
    if log_output then
        [fd, err] = mopen("log_adams_bashforth_passo_4.txt", "a+")
        if err ~= 0 then
            return
        end
    end
    
    if log_output then
        mfprintf(fd, "%f,%f,%f,%f,%f,%f,%s\n", a,b,N,h,x0,y0,method)
        mfprintf(fd, "x_i,y_i_aproximado,y_i_analitico,erro\n")
    end
    
    tic()
    
    disp("METODO DE PASSO UNICO")
    
    [x_aux, y_aux, tempo_aux] = metodo_passo_unico(a, b, 3, h, x0, y0, method)
    
    disp("METODO DE PASSO MULTIPLO")
    
    x = zeros(N+1)
    y = zeros(N+1)
    y_analitico = zeros(N+1)
    x(1:4) = x_aux(1:4)
    y(1:4) = y_aux(1:4)
    y_analitico(1:4) = f_solucao(x_aux(1:4))  
    
    if log_output then
        for i = 2:4
            erro = abs(y_analitico(i) - y(i))
            mfprintf(fd, "%.16f,%.16f,%.16f,%.16f\n", x(i), y(i), y_analitico(i), erro)
        end
    end 
    
    
    for i = 4:N
        
        x(i+1) = x(i) + h
        y(i+1) = y(i) + (h/24)*(55*f_pvi(x(i), y(i)) - 59*f_pvi(x(i-1), y(i-1)) + 37*f_pvi(x(i-2), y(i-2)) - 9*f_pvi(x(i-3), y(i-3)))
        y_analitico(i+1) = f_solucao(x(i+1))
        erro = abs(y_analitico(i+1) - y(i+1))
        
        disp([x(i+1) y(i+1) y_analitico(i+1) erro])
        
        if log_output then
            mfprintf(fd, "%.16f,%.16f,%.16f,%.16f\n", x(i+1), y(i+1), y_analitico(i+1), erro)
        end
        
    end
    
    tempo = toc()
    
    if log_output then
        mfprintf(fd, "Tempo=%f\n\n", tempo)
        mclose(fd);
    end
    
    clf
    plot2d(x', [y y_analitico], style=[20 1], leg="PVI@ANALITICO")
    title("SOLUCAO", "color", "blue")
    xlabel("x", "color", "blue")
    ylabel("y", "color", "blue")
    
endfunction



function [x, y, tempo] = metodo_adams_bashforth(p, a, b, N, h, x0, y0, method, log_output)
    
    if ~exists('log_output', 'local') then
        log_output = %T
    end
    
    select p
    case 2 then
        [x, y, tempo] = metodo_adams_bashforth_passo_2(a, b, N, h, x0, y0, method, log_output)
    case 3 then
        [x, y, tempo] = metodo_adams_bashforth_passo_3(a, b, N, h, x0, y0, method, log_output)
    case 4 then
        [x, y, tempo] = metodo_adams_bashforth_passo_4(a, b, N, h, x0, y0, method, log_output)
    else
        break
    end
    
endfunction



// LaTeX $x_{n+1} = x_{n} + h$
// LaTeX $y_{n+1} = y_n + \frac{h}{2} \cdot (f(x_{n+1}, y_{n+1}) + f(x_n, y_n))$
function [x, y, tempo] = metodo_adams_moulton_passo_2(a, b, N, h, x0, y0, method, log_output)
    
    disp("METODO DE ADAMS MOULTON DE PASSO 2")
    
    if ~exists('log_output', 'local') then
        log_output = %T
    end
    
    if log_output then
        [fd, err] = mopen("log_adams_moulton_passo_2.txt", "a+")
        if err ~= 0 then
            return
        end
    end
    
    if log_output then
        mfprintf(fd, "%f,%f,%f,%f,%f,%f,%s\n", a,b,N,h,x0,y0,method)
        mfprintf(fd, "x_i,y_i_aproximado,y_i_analitico,erro\n")
    end
    
    tic()
    
    disp("METODO DE PASSO UNICO")
    
    [x_aux, y_aux, tempo_aux] = metodo_passo_unico(a, b, N, h, x0, y0, method)
    
    disp("METODO DE PASSO MULTIPLO")
    
    x = zeros(N+1)
    y = zeros(N+1)
    y_analitico = zeros(N+1)
    x(1) = x0
    y(1) = y0
    y_analitico(1) = y0
    
    for i = 1:N
        
        x(i+1) = x(i) + h
        y(i+1) = y(i) + (h/2)*(f_pvi(x(i+1), y_aux(i+1)) + f_pvi(x(i), y(i)))
        y_analitico(i+1) = f_solucao(x(i+1))
        erro = abs(y_analitico(i+1) - y(i+1))
        
        disp([x(i+1) y(i+1) y_analitico(i+1) erro])
        if log_output then
            mfprintf(fd, "%.16f,%.16f,%.16f,%.16f\n", x(i+1), y(i+1), y_analitico(i+1), erro)
        end
        
    end
    
    tempo = toc()
    
    if log_output then
        mfprintf(fd, "Tempo=%f\n\n", tempo)
        mclose(fd);
    end
    
    clf
    plot2d(x', [y y_analitico], style=[20 1], leg="PVI@ANALITICO")
    title("SOLUCAO", "color", "blue")
    xlabel("x", "color", "blue")
    ylabel("y", "color", "blue")
    
endfunction



// LaTeX $x_{n+1} = x_{n} + h$
// LaTeX $y_{n+1} = y_n + \frac{h}{12} \cdot (5 \cdot f(x_{n+1}, y_{n+1}) + 8 \cdot f(x_n, y_n) - f(x_{n-1}, y_{n-1}))$
function [x, y, tempo] = metodo_adams_moulton_passo_3(a, b, N, h, x0, y0, method, log_output)
    
    disp("METODO DE ADAMS MOULTON DE PASSO 3")
    
    if ~exists('log_output', 'local') then
        log_output = %T
    end
    
    if log_output then
        [fd, err] = mopen("log_adams_moulton_passo_3.txt", "a+")
        if err ~= 0 then
            return
        end
    end
    
    if log_output then
        mfprintf(fd, "%f,%f,%f,%f,%f,%f,%s\n", a,b,N,h,x0,y0,method)
        mfprintf(fd, "x_i,y_i_aproximado,y_i_analitico,erro\n")
    end
    
    tic()
    
    disp("METODO DE PASSO UNICO")
    
    [x_aux, y_aux, tempo_aux] = metodo_passo_unico(a, b, N, h, x0, y0, method)
    
    disp("METODO DE PASSO MULTIPLO")
    
    x = zeros(N+1)
    y = zeros(N+1)
    y_analitico = zeros(N+1)
    x(1:2) = x_aux(1:2)
    y(1:2) = y_aux(1:2)
    y_analitico = f_solucao(x_aux(1:2))
    
    if log_output then
        for i = 2:2
            erro = abs(y_analitico(i) - y(i))
            mfprintf(fd, "%.16f,%.16f,%.16f,%.16f\n", x(i), y(i), y_analitico(i), erro)
        end
    end 
    
    for i = 2:N
        
        x(i+1) = x(i) + h
        y(i+1) = y(i) + (h/12)*(5*f_pvi(x(i+1), y_aux(i+1)) + 8*f_pvi(x(i), y(i)) - f_pvi(x(i-1), y(i-1)))
        y_analitico(i+1) = f_solucao(x(i+1))
        erro = abs(y_analitico(i+1) - y(i+1))
        
        disp([x(i+1) y(i+1) y_analitico(i+1) erro])
        if log_output then
            mfprintf(fd, "%.16f,%.16f,%.16f,%.16f\n", x(i+1), y(i+1), y_analitico(i+1), erro)
        end
        
    end
    
    tempo = toc()
    
    if log_output then
        mfprintf(fd, "Tempo=%f\n\n", tempo)
        mclose(fd);
    end
    
    clf
    plot2d(x', [y y_analitico], style=[20 1], leg="PVI@ANALITICO")
    title("SOLUCAO", "color", "blue")
    xlabel("x", "color", "blue")
    ylabel("y", "color", "blue")
    
endfunction



// LaTeX $x_{n+1} = x_{n} + h$
// LaTeX $y_{n+1} = y_{n} + \frac{h}{24} \cdot (9 \cdot f(x_{n+1}, y_{n+1}) + 19 \cdot f(x_n, y_n) - 5 \cdot f(x_{n-1}, y_{n-1}) + f(x_{n-2}, y_{n-2}))$
function [x, y, tempo] = metodo_adams_moulton_passo_4(a, b, N, h, x0, y0, method, log_output)
    
    disp("METODO DE ADAMS MOULTON DE PASSO 4")
    
    if ~exists('log_output', 'local') then
        log_output = %T
    end
    
    if log_output then
        [fd, err] = mopen("log_adams_moulton_passo_4.txt", "a+")
        if err ~= 0 then
            return
        end
    end
    
    if log_output then
        mfprintf(fd, "%f,%f,%f,%f,%f,%f,%s\n", a,b,N,h,x0,y0,method)
        mfprintf(fd, "x_i,y_i_aproximado,y_i_analitico,erro\n")
    end
    
    tic()
    
    disp("METODO DE PASSO UNICO")
    
    [x_aux, y_aux, tempo_aux] = metodo_passo_unico(a, b, N, h, x0, y0, method)
    
    disp("METODO DE PASSO MULTIPLO")
    
    x = zeros(N+1)
    y = zeros(N+1)
    y_analitico = zeros(N+1)
    x(1:3) = x_aux(1:3)
    y(1:3) = y_aux(1:3)
    y_analitico(1:3) = f_solucao(x_aux(1:3))
    
    if log_output then
        for i = 2:3
            erro = abs(y_analitico(i) - y(i))
            mfprintf(fd, "%.16f,%.16f,%.16f,%.16f\n", x(i), y(i), y_analitico(i), erro)
        end
    end 
    
    for i = 3:N
        
        x(i+1) = x(i) + h
        y(i+1) = y(i) + (h/24)*(9*f_pvi(x(i+1), y_aux(i+1)) + 19*f_pvi(x(i), y(i)) - 5*f_pvi(x(i-1), y(i-1)) + f_pvi(x(i-2), y(i-2)))
        y_analitico(i+1) = f_solucao(x(i+1))
        erro = abs(y_analitico(i+1) - y(i+1))
        
        disp([x(i+1) y(i+1) y_analitico(i+1) erro])
        if log_output then
            mfprintf(fd, "%.16f,%.16f,%.16f,%.16f\n", x(i+1), y(i+1), y_analitico(i+1), erro)
        end
        
    end
    
    tempo = toc()
    
    if log_output then
        mfprintf(fd, "Tempo=%f\n\n", tempo)
        mclose(fd);
    end
    
    clf
    plot2d(x', [y y_analitico], style=[20 1], leg="PVI@ANALITICO")
    title("SOLUCAO", "color", "blue")
    xlabel("x", "color", "blue")
    ylabel("y", "color", "blue")
    
endfunction



function [x, y, tempo] = metodo_adams_moulton(p, a, b, N, h, x0, y0, method, log_output)
    
    if ~exists('log_output', 'local') then
        log_output = %T
    end
    
    select p
    case 2 then
        [x, y, tempo] = metodo_adams_moulton_passo_2(a, b, N, h, x0, y0, method, log_output)
    case 3 then
        [x, y, tempo] = metodo_adams_moulton_passo_3(a, b, N, h, x0, y0, method, log_output)
    case 4 then
        [x, y, tempo] = metodo_adams_moulton_passo_4(a, b, N, h, x0, y0, method, log_output)
    else
        break
    end
    
endfunction



function [x,y,tempo] = metodo_preditor_corretor(a, b, N, h, x0, y0, tol, log_output)
    
    function [r] = f_pred(i)
        r = y(i) + h*f_pvi(x(i),y(i))
    endfunction
    
    function [r] = f_corr(i, pred)
                
        while 1
            
            y_new = y(i) + (h/2)*(f_pvi(x(i), y(i)) + f_pvi(x(i+1),pred))
            erro = abs(y_new - pred)/abs(y_new)
            
            if erro < tol
                r = y_new
                break
            end
            
            pred = y_new
            
        end
        
    endfunction
    
    if ~exists('log_output', 'local') then
        log_output = %T
    end
    
    if log_output then
        [fd, err] = mopen("log_preditor_corretor.txt", "a+")
        if err ~= 0 then
            return
        end
    end
    
    if log_output then
        mfprintf(fd, "%f,%f,%f,%f,%f,%f,%f\n", a,b,N,h,x0,y0,tol)
        mfprintf(fd, "x_i,y_i_aproximado,y_i_analitico,erro\n")
    end
    
    x = zeros(N+1)
    y = zeros(N+1)
    y_analitico = zeros(N+1)
    x(1) = x0
    y(1) = y0
    y_analitico(1) = y0
    
    tic()
    
    disp("METODO PREDITOR CORRETOR")    
    
    for i = 1:N
        
        x(i+1) = x(i) + h
        pred_i = f_pred(i)
        corr_i = f_corr(i, pred_i)
        y(i+1) = corr_i
        y_analitico(i+1) = f_solucao(x(i+1))
        erro = abs(y_analitico(i+1) - y(i+1))
        
        disp([x(i+1) y(i+1) y_analitico(i+1) erro])
        if log_output then
            mfprintf(fd, "%.16f,%.16f,%.16f,%.16f\n", x(i+1), y(i+1), y_analitico(i+1), erro)
        end
        
    end
    
    tempo = toc()
    

    if log_output then
        mfprintf(fd, "Tempo=%f\n\n", tempo)
        mclose(fd);
    end
    
    clf
    plot2d(x', [y y_analitico], style=[20 1], leg="PVI@ANALITICO")
    title("SOLUCAO", "color", "blue")
    xlabel("x", "color", "blue")
    ylabel("y", "color", "blue")
    
endfunction



function [x,y,tempo] = metodo_preditor_corretor_implicito(a, b, N, h, x0, y0, tol, method, N_explicito, log_output)
    
    function [r] = f_pred(i)
        r = y(i-1) + h*f_pvi(x(i-1),y(i-1))
    endfunction
    
    function [r] = f_corr(i, pred)
                
        while 1
            
            y_new = y(i-1) + (h/3)*(f_pvi(x(i-1), y(i-1)) + 4 * f_pvi(x(i), y(i)) + f_pvi(x(i+1), pred))
            erro = abs(y_new - pred)/abs(y_new)
            
            if erro < tol
                r = y_new
                break
            end
            
            pred = y_new
            
        end
        
    endfunction
    
    if ~exists('log_output', 'local') then
        log_output = %T
    end
    
    if log_output then
        [fd, err] = mopen("log_preditor_corretor_implicito.txt", "a+")
        if err ~= 0 then
            return
        end
    end
    
    if log_output then
        mfprintf(fd, "%f,%f,%f,%f,%f,%f,%f,%s\n", a,b,N,h,x0,y0,tol,method)
        mfprintf(fd, "x_i,y_i_aproximado,y_i_analitico,erro\n")
    end
    
    tic()
    
    disp("METODO DE PASSO UNICO")
    
    [x_aux, y_aux, tempo_aux] = metodo_passo_unico(a, b, N_explicito, h, x0, y0, method)
    
    disp("METODO PREDITOR CORRETOR")
    
    x = zeros(N+1)
    y = zeros(N+1)
    y_analitico = zeros(N+1)
    x(1:N_explicito+1) = x_aux(1:N_explicito+1)
    y(1:N_explicito+1) = y_aux(1:N_explicito+1)
    y_analitico(1:N_explicito+1) = f_solucao(x_aux(1:N_explicito+1))
    
    if log_output then
        for i = 2:(N_explicito+1)
            erro = abs(y_analitico(i) - y(i))
            mfprintf(fd, "%.16f,%.16f,%.16f,%.16f\n", x(i), y(i), y_analitico(i), erro)
        end
    end 
    
    for i = N_explicito+1:N
        
        x(i+1) = x(i) + h
        pred_i = f_pred(i)
        corr_i = f_corr(i, pred_i)
        y(i+1) = corr_i
        y_analitico(i+1) = f_solucao(x(i+1))
        erro = abs(y_analitico(i+1) - y(i+1))
        
        disp([x(i+1) y(i+1) y_analitico(i+1) erro])
        if log_output then
            mfprintf(fd, "%.16f,%.16f,%.16f,%.16f\n", x(i+1), y(i+1), y_analitico(i+1), erro)
        end
        
    end
    
    tempo = toc()
    
    if log_output then
        mfprintf(fd, "Tempo=%f\n\n", tempo)
        mclose(fd);
    end
    
    clf
    plot2d(x', [y y_analitico], style=[20 1], leg="PVI@ANALITICO")
    title("SOLUCAO", "color", "blue")
    xlabel("x", "color", "blue")
    ylabel("y", "color", "blue")
    
    
endfunction






















