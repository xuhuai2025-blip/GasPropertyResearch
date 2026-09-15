function value=sha256(path)
%SHA256 Fingerprint installed source files without copying licensed data.
fid=fopen(path,'rb');
if fid<0,error('refpropbuild:MissingFile','Cannot read %s.',path);end
cleanup=onCleanup(@()fclose(fid));
bytes=fread(fid,Inf,'*uint8');
digest=java.security.MessageDigest.getInstance('SHA-256');
digest.update(bytes);
value=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end
