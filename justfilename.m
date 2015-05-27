function fn = justfilename(filenames,remove_extension)

if(nargin < 2)
    remove_extension = false;
end

is1 = ~iscell(filenames);
if(~iscell(filenames))
    filenames = {filenames};
end
fn = filenames;

if(ispc)
    hasdir = cellfun(@(x)(any(x=='/' | x=='\')),filenames);

    fn(~hasdir) = filenames(~hasdir);
    fn(hasdir) = cellfun(@(x)(x(max(find(x=='/' | x=='\'))+1:end)),filenames(hasdir),'uniformoutput',false);
else
    hasdir = cellfun(@(x)(any(x=='/')),filenames);

    fn(~hasdir) = filenames(~hasdir);
    fn(hasdir) = cellfun(@(x)(x(max(find(x=='/'))+1:end)),filenames(hasdir),'uniformoutput',false);    
end

if(remove_extension)
    hasext = cellfun(@(x)(any(x=='.')),fn);
    fn(hasext) = cellfun(@(x)(x(1:max(find(x=='.'))-1)),fn(hasext),'uniformoutput',false);
end

if(is1)
    fn = fn{1};
end