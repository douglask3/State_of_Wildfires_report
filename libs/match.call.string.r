match.call.string <- function(n=1,call=sys.call(sys.parent(n = n)),...) {
	FunName = as.character(match.call(call=call,...))[1]
	if (FunName == "match.call.string")  FunName = "N/A"
	return(FunName)
}
	