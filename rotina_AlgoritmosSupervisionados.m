y=getGrupos('iris.txt');

n_grupos=length(y(end,:));


cols = length(y(1).mat(end,:));
 
 
 % Vai fazer um loop percorrendo todos os grupos;
for ngrp=1:n_grupos
    % O n_esimo grupo escolhido para fazer os teste com o algoritmo
    % supervisionado;
    grupo=y(ngrp).mat;
    for cols_mat=1: (length(grupo(end,:)))
        % Aqui é selecionado a coluna que será classe, inicialmente a
        % coluna 1 ate a ultima coluna da matriz;
        classe=cols_mat;
        % Funcao retorna nesta variavel as colunas que não sao classes da
        % matria;
        vet_col_semClasse=getMatrizSemClasseX(grupo,classe);
        
        mdl=fitcnb(grupo(:,vet_col_semClasse),grupo(:,classe),'DistributionNames','mn');
        
        
        
       
    end
    
    
   
end