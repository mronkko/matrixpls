#'@return
#'
<%

attributeFile <- read.csv2("man-roxygen/attributes.csv",stringsAsFactors = FALSE)
attributes <- unlist(strsplit(attributes,","))

%>
#'\code{<%=attributes[1]%>} returns the following as attributes:
#'
<%
for(attribute in unlist(strsplit(attributes[2]," "))){
%>
#'\item{<%=attribute%>}{<%=attributeFile$description[attributeFile$attribute == attribute]%>.}
<%
}

if(length(attributes) == 3){
%>
#'@return
#'Additionally, all attributes returned by functions called by <%=attributes[1]%> are returned.
#'This can include:
<%
  for(attribute in unlist(strsplit(attributes[3]," "))){
%>
#'\item{<%=attribute%>}{<%=attributeFile$description[attributeFile$attribute == attribute]%>.}
<%
  }
}
%>
