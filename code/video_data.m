function [V, E, fid, ground_truth] = video_data(vname,phi,E, fid_perc, train_test_split, fid_sampling, varargin)

V = phi;
E = E;

%[N, Neig, NClass] = get_dims(vname,V);
if (numel(varargin) == 1)
    
    [fid, ground_truth] = get_fidelity(vname, fid_perc, train_test_split, fid_sampling, varargin{1});
else
    [fid, ground_truth] = get_fidelity(vname, fid_perc, train_test_split,fid_sampling);
end


end


function [fid, gt] = get_fidelity(vname, fid_perc, train_test_split, fid_sampling, varargin)



    gt = load(['data/gt_folder/',vname,'.mat']);
    gt = gt.gt;
    display(size(gt))
    T_segment = 6;
    Nclass = 15;
    
    
    fid = {};
    
    % Choose 2 Slices for each class
    
    N = size(gt,1);
    
    
    
    desired_label = [1,2,3,4,5,6,7,9,10,11,12,13,14,15,16];
    
    
    % Training and testing split
    cutoff = round( size(gt,1) * train_test_split)
    temp = gt(1:cutoff);
       
    
    %% Fidelity Sampling Method
    
     % Sample fidelity within each category
     if (fid_sampling == 1 || fid_sampling == 2)
             
         for i = 1:Nclass
 
                L = desired_label(i);
 
                vec = find(temp == L);
                %vec = vec(1:round(numel(vec)/2));
 
                % Even sample
                n = length(vec);
 
                if n == 0
                    fid{L} = [];
                else
                    number = round(n * fid_perc);
                    if fid_sampling == 1 % Uniform within each category
                        fid{L} = vec(round(linspace(1,n,number)));
                    elseif fid_sampling == 2 % Random within each category
                        index = randi(n,[number,1]);
                        fid{L} = vec(index);
                    end
                    
                end
 
        end
     end
    
        
    if (fid_sampling == 3)
        index = randi(numel(temp),[round(fid_perc * numel(temp)),1]);
        t = temp(index);
                 
        for i = 1:Nclass

           L = desired_label(i);

           vec = find(t == L);
           n = length(vec);
           if n == 0
               v = find(temp==L);
               
               if numel(v)>0
                   fid{L} = [v(1)];
               end
               
           else
               fid{L} = vec;
   
           end

        end
    end
        
        
end
    
    

    
    
    
    
    
    



