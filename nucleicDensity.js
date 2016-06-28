
urlstr=document.location.pathname;
var r = /[^\/]*$/;
urlstr=urlstr.replace(r, ''); 
urlstr=urlstr.replace(/^\//, '');

function outputDensityUI(){
	output='<div style="font-style:bold;margin:0px;margin-top:30px;margin-left:10px;background-color:#f0f0f0;"> \
        Additional Tools:\
    </div> \
    <div style="margin-left:10px;border:1px solid #606060;"> \
    <input id="btnCallWebServiceDensity" value="Compute A, T, G, C, G+C, A+T" type="button" />  \
    <input id="btnCallGCSkew" value="Generate GC Skew" type="button" />  \
	<div style="margin-bottom:10px;"> \
    start: <input id="sbegin" type="text" /> press \'q\' on keyboard to set<br /> \
    end: <input id="send" type="text" /> press \'z\' on keyboard to set<br /> \
	length: <span id="sequenceLength"></span> nucleotides\
	</div> \
	    <div id="status"></div> \
    <div id="result"> \
    </div> <br clear="left"> \
    <span style="margin-left:10px;background-color:#f0f0f0;">Result Log:</span> \
    <div id="outfile"></div> \
  	</div> \
	';
	
	document.write(output);
}

function otherCredits(){
		//no other credits
}
 
        $(document).ready(function () {
									
				$('#sbegin').change(function() {
  						updateSequenceLength();
				});			
				$('#send').change(function() {
  						updateSequenceLength();
				});	
									
        	 	$("#refseq").val(usa);
                $("#sbegin").val(sbegin);
                $("#send").val(send);
				
				sequenceLength=parseInt(send)-parseInt(sbegin)+1;
				$("#sequenceLength").html(sequenceLength);
                
                
                $('body').keyup(function (event) {
						    var direction = null;
						    // handle keys
						    if (event.keyCode == 81) {
						      // set start
						      direction=$("#Nucleotide").html();
						      $("#sbegin").val(direction);
							  updateSequenceLength();
						    } else if (event.keyCode == 90) {
						      // set end
						      direction=$("#Nucleotide").html();
						      $("#send").val(direction);
							  updateSequenceLength();
						    }
						  	});
                
                
                
			    $("#btnCallWebServiceDensity").click(function (event) {
			    		
            	 	
            		sbegin=$("#sbegin").val();
                send=$("#send").val();
                if (parseInt(sbegin)>=parseInt(send)) {
            	  	alert('invalid input: start of sequence should be smaller than end');
            	  	return;
            	  	}
            		
            		sbegin=$("#sbegin").val();
            		send=$("#send").val();
								sequenceLength=$("#sequenceLength").html();
                var wsUrl = "../../density.php";
				
								uiLoading ("Computing...");
				
                $.ajax({
                    type: "GET",
                    data: {
         									  start: sbegin,
            								end: send,
            								filepath: urlstr
        								},
                    url: wsUrl,
                    contentType: "text/html",
                    success: processSuccessDensity,
                    error: processError
                });
								
            
          
        });
		});
		
	
		function updateSequenceLength(){
			sbegin=$("#sbegin").val();
      send=$("#send").val();
			sequenceLength=parseInt(send)-parseInt(sbegin)+1;
			$("#sequenceLength").html(sequenceLength);
		}
        
        
    function uiLoading (message) {
    	//var bufferOutFile=$("#outfile").val();
    	$("#status").html(" <img src='../../loading.gif' style='float:left;' /><div style='font-size:11pt;padding-top:10px;padding-bottom:15px;'>"+message+"</div>" );
    	
    }

	
		 function processSuccessDensity(response) {
						  $("#outfile").prepend("<div class='densityTable'>"+response+"</div>");
				      $("#status").html("Completed." );
        }

        function processError() {
        	
            $("#outfile").prepend("<br />Error.  Connection problem. <div class='resultdivider'></div>");
            $("#status").html("Completed with error." );
        }  

