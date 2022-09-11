
% TP3 de Probabilites : fonctions a completer et rendre sur Moodle
% Nom : Foucher
% Prénom : Nathan
% Groupe : 1SN-C

function varargout = fonctions_TP3_proba(varargin)

    switch varargin{1}
        
        case 'matrice_inertie'
            
            [varargout{1},varargout{2}] = feval(varargin{1},varargin{2:end});
            
        case {'ensemble_E_recursif','calcul_proba'}
            
            [varargout{1},varargout{2},varargout{3}] = feval(varargin{1},varargin{2:end});
    
    end
end

% Fonction ensemble_E_recursif (exercie_1.m) ------------------------------
function [E,contour,G_somme] = ...
    ensemble_E_recursif(E,contour,G_somme,i,j,voisins,G_x,G_y,card_max,cos_alpha)

    contour(i,j)=0;
    for k=1:size(voisins) 
        i_voisin=voisins(k,1)+i;
        j_voisin=voisins(k,2)+j;
        G_voisin = [G_x(i_voisin,j_voisin),G_y(i_voisin, j_voisin)];
        
        if contour(i_voisin,j_voisin) ~= 0 && dot(G_voisin/norm(G_voisin),(G_somme/norm(G_somme))) >= cos_alpha
            G_somme(1) = G_somme(1) + G_x(i_voisin,j_voisin);
            G_somme(2) = G_somme(2) + G_y(i_voisin,j_voisin);
            E=[E; i_voisin,j_voisin];
            [E,contour,G_somme]=ensemble_E_recursif(E,contour,G_somme,i_voisin,j_voisin,voisins,G_x,G_y,card_max,cos_alpha);
             end

        end 
 
end

% Fonction matrice_inertie (exercice_2.m) ---------------------------------
function [M_inertie,C] = matrice_inertie(E,G_norme_E)
    n=length(E);
    G_norme_new=G_norme_E/sum(G_norme_E);
    C=(G_norme_new'*fliplr(E));
    M11=(1/sum(G_norme_E))*((E(:,2)-C(1).*ones(n,1)).^2)'*G_norme_E;
    M12=(1/sum(G_norme_E))*((E(:,2)-C(1).*ones(n,1)).*(E(:,1)-C(2)*ones(n,1)))'*G_norme_E;
    M22=(1/sum(G_norme_E))*((E(:,1)-C(2)*ones(n,1)).^2)'*G_norme_E;
    M_inertie=[M11, M12 ;M12 , M22];
end

% Fonction calcul_proba (exercice_2.m) ------------------------------------
function [x_min,x_max,probabilite] = calcul_proba(E_nouveau_repere,p)
    x_min=min(E_nouveau_repere(:,1));
    y_min=min(E_nouveau_repere(:,2));
    x_max=max(E_nouveau_repere(:,1));
    y_max=max(E_nouveau_repere(:,2));
    
    N=round((y_max-y_min)*(x_max-x_min));
    n=length(E_nouveau_repere);
    probabilite=1-binocdf(n-1,N,p);
    
end
