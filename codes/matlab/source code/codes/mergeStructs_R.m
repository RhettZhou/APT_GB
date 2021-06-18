function [merged_struct] = mergeStructs_R(struct_a,struct_b)
%%if one of the structres is empty do not merge
if isempty(struct_a)
    merged_struct=struct_b;
    return
end
if isempty(struct_b)
    merged_struct=struct_a;
    return
end
%%insert struct a
merged_struct=struct_a;
f = fieldnames(merged_struct);
%%insert struct b
Vsize_A = length(struct_a(1).(f{1}));
for i = 1:length(f)
    size_a = length(struct_a(1).(f{i}));
    size_b = length(struct_b(1).(f{i}));
    for j=1:size_b
        if i == 2
            merged_struct(1).(f{i})(size_a+j,:) = struct_b.(f{i})(j,:) + Vsize_A*ones(1,3);
        else
            merged_struct(1).(f{i})(size_a+j,:) = struct_b.(f{i})(j,:);
        end
    end
end