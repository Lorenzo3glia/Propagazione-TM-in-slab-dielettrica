function[zero]=findzero(funfcn,a1,b1,acc,varargin);

zero=[];
fa=feval(funfcn,a1,varargin{:});
fb=feval(funfcn,b1,varargin{:});
if (fa==0), zero=a1;
elseif (fb==0), zero=b1;
elseif (fa*fb<0),
   c1=(a1+b1)/2;
   fc=feval(funfcn,c1,varargin{:});
   while(abs(fc)>acc),
      c1=(a1+b1)/2;
      fc=feval(funfcn,c1,varargin{:});
      if (fa*fc<0),
         b1=c1;
      else
         a1=c1;
      end
   end
   zero=c1;
end


%La funzione feval restituisce il valore di della funzione in un punto insrito come dato

%La funzione non fa altro che implementare un semplice metodo di bisezione,
%nello specifico, in input viene fornita la funzione di cui si vogliono
%valutare gli zeri e gli estremi di valutazione con la relativa
%accuratezza.Vengono inseriti degli eventuali paramentri aggiuntivi
%attraverso il parametro varargin.
%Nelle parti successive non si fa altro che implementare un metodo di
%bisezione con cui valutare gli zeri della funzione.

