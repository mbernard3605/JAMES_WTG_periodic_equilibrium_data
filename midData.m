function [ smallerData ] = midData( fullData )

smallerData=(fullData(2:end,:)+fullData(1:end-1,:))/2;

end

