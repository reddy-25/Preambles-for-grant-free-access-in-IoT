% Generating premables using ZadoffChu Sequence

function [preambles_set ,para] = gen_preambs_common(para)

preambles_set = [];
zadoff_roots = para.no_of_pre;

if para.long_preamble == 1
    for ii = 1:zadoff_roots      % 1:No of roots
        preambles_set(:,ii) = zadoffChuSeq(ii,para.pre_length);
    end
    
elseif para.short_preamble== 1
    num_zeros = para.num_preamble_sc - para.short_pre_length;
    for ii = 1:zadoff_roots      % 1:No of roots
        preambles_set_1 = [zadoffChuSeq(ii,para.short_pre_length); zeros(num_zeros,1)];
        preambles_set_1 = repmat(preambles_set_1,1,para.time_sym);
        preambles_set(:,ii) = preambles_set_1(:);
    end
elseif para.spread_sequence==1
    %for multiple preambles case with spread sequence
    num_zeros = para.num_preamble_sc - para.short_pre_length;
    %code=hadamard(para.time_sym);
    code = zadoffChuSeq(1,7);
    
    code1 = code;
    
    temp=1;
    
    %generating 2 preambles 
    for ii=1:2
        
          preambles_set_1 = [zadoffChuSeq(ii,para.short_pre_length); zeros(num_zeros,1)];
          
          preambles_set_1 = repmat(preambles_set_1,1,para.time_sym);
          preambles_set_3(:,ii) = preambles_set_1(:);
          preambles_set_2 = preambles_set_1;
        
          
          for ll=1:6
              code = circshift(code1,ll-1);
              preambles_set_1=preambles_set_2;
              %multiplying symbol wise with each code sequence
              for jj= 1:para.time_sym
                 
                  preambles_set_1(:,jj) = code(jj,1).*preambles_set_1(:,jj);
                  
              end
              preambles_set(:,temp) = preambles_set_1(:);
              temp=temp+1;
          end
    end
end

end