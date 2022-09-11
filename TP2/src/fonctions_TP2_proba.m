
% TP2 de Probabilites : fonctions a completer et rendre sur Moodle
% Nom : Foucher 
% Prenom : Nathan
% Groupe : 1SN-C

function varargout = fonctions_TP2_proba(varargin)

    switch varargin{1}
        
        case {'calcul_frequences_caracteres','determination_langue','coeff_compression','gain_compression','partitionnement_frequences'}

            varargout{1} = feval(varargin{1},varargin{2:end});
            
        case {'recuperation_caracteres_presents','tri_decroissant_frequences','codage_arithmetique'}
            
            [varargout{1},varargout{2}] = feval(varargin{1},varargin{2:end});
            
        case 'calcul_parametres_correlation'
            
            [varargout{1},varargout{2},varargout{3}] = feval(varargin{1},varargin{2:end});
            
    end

end

% Fonction calcul_frequences_caracteres (exercice_0.m) --------------------
function frequences = calcul_frequences_caracteres(texte,alphabet)
    frequences=zeros(length(alphabet),1);
    for i=1:length(alphabet)
        frequences(i)=length(find(texte==alphabet(i)));
    end
    frequences=frequences/sum(frequences);
 end

% Fonction recuperation_caracteres_presents (exercice_0.m) ----------------
function [selection_frequences,selection_alphabet] = ...
                      recuperation_caracteres_presents(frequences,alphabet)

    indices=find(frequences>0);
    selection_alphabet=alphabet(indices);
    selection_frequences=frequences(indices);

end

% Fonction tri_decremental_frequences (exercice_0.m) ----------------------
function [frequences_triees,indices_frequences_triees] = ...
                           tri_decroissant_frequences(selection_frequences)
    [frequences_triees,indices_frequences_triees] = sort(selection_frequences,'descend');
end

% Fonction determination_langue (exercice_1.m) ----------------------------
function langue = determination_langue(frequences_texte, frequences_langue, nom_norme)
    % Note : la variable nom_norme peut valoir 'L1', 'L2' ou 'Linf'.

    tfrequences_texte=transpose(frequences_texte);
    mat_freq=repmat(tfrequences_texte,15,1);

    switch nom_norme
        case 'L2'
            mat_error=(mat_freq-frequences_langue).^2;
            error=sum(mat_error,2);
        case 'L1'
            mat_error=abs(mat_freq-frequences_langue);
            error=sum(mat_error,2);
        case 'Linf'
            mat_error=abs(mat_freq-frequences_langue);
            error=max(mat_error,[],2);
    end

    [x,y]=find(error==min(error));
    langue=x;
end

% Fonction coeff_compression (exercice_2.m) -------------------------------
function coeff_comp = coeff_compression(signal_non_encode,signal_encode)
    coeff_comp=size(signal_non_encode,1)*size(signal_non_encode,2)*8/(size(signal_encode,1)*size(signal_encode,2));
end

% Fonction coeff_compression (exercice_2bis.m) -------------------------------
function gain_comp = gain_compression(coeff_comp_avant,coeff_comp_apres)
    gain_comp=coeff_comp_apres/coeff_comp_avant;
end

% Fonction partitionnement_frequences (exercice_3.m) ----------------------
function bornes = partitionnement_frequences(selection_frequences)



end

% Fonction codage_arithmetique (exercice_3.m) -----------------------------
function [borne_inf,borne_sup] = ...
                       codage_arithmetique(texte,selection_alphabet,bornes)


    
end