%% preGRPDA-L
ga_pre_list=[1/10,1/5,1,5,10]; 
line_clr=["r","g","b","r--","g--","b--"];
iter=size(ga_pre_list,2);
f_pregrpda_L=cell(1,iter);
t_pregrpda_L=cell(1,iter);
fstar=inf;

%solver
for i=1:iter
    a=(2-sqrt(2))/2;  
    ga=ga_pre_list(i);
    h_a=(2*a*a-6*a+3)/((1-a)*(a*a-3*a+2));
    step=h_a;   
    gama1=sqrt(1/step)/ga;
    gama2=sqrt(1/step)*ga;
    M1=gama1*sum(abs(K))';
    M2A=gama2*sum(abs(A),2);
    M2B=gama2*sum(abs(B),2);
    
    a=0.5; 
    kesai=1.1;
    [x_pregrpda_L,f_pregrpda_L{i},t_pregrpda_L{i}]= preGRPDA_L(A,B,b,a,M1,M2A,M2B,kesai);

    fstar= min(fstar,min(f_pregrpda_L{i}));
end

%plot
for i=1:iter
    figure(1)
    semilogy(t_pregrpda_L{i},(f_pregrpda_L{i}-optval)/optval,line_clr(i),'DisplayName',num2str(ga_pre_list(i)));
    hold on
    legend show;
end
