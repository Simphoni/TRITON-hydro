$('#configGenForm').on("submit", function(event) {
    //stop the form from fully submitting and refreshing the page
	event.preventDefault();
	//get the file from input
	//var file = $("#XMLFilename").prop('files');
	var file = handleFileSelect('XMLFilename');
	getFileContent(file);
	//content = file.content();
	//validate the file
	//TO DO
	//parse the xml
	//asyncFunc();
	//console.log(content);
})

function handleFileSelect(fileInput) {
	if (!window.File || !window.FileReader || !window.FileList || !window.Blob) {
		alert('The File APIs are not fully supported in this browser.');
		return;
	}   

	var input = document.getElementById(fileInput);
	if (!input) {
		alert("Um, couldn't find the fileinput element.");
	}
	else if (!input.files) {
		alert("This browser doesn't seem to support the `files` property of file inputs.");
	}
	else if (!input.files[0]) {
		alert("Please select a file before clicking 'Generate Form'");               
	}
	else {
		var file = input.files[0];
		//get the filename and put it into cover text area
		return file;
		//var fr = new FileReader();
		//fr.onload = receivedText;
		//var text = fr.readAsText(file);
		//return text;
		//fr.readAsDataURL(file);
	}
}

async function getFileContent(file, onLoadCallback) {
	var fr = new FileReader();
	var text;
	fr.onload = function() {
		//I have come to the conclusion
		//it is impossible to get the value
		//of text out of this block, so you
		//have to do ALL the string handling
		//right here, neat
		text = fr.result;
		handleXML(text);
	};
	fr.readAsText(file);
}

function handleXML(str) {
	//console.log(str);
	$xml = $($.parseXML(str));
	//get the root node, should always be 'form'
	$root = $xml.find(":root");
	var form = document.createElement('form');
	form.setAttribute("id", $($root).attr('id'));
	form.setAttribute("method", $($root).attr('method'));
	form.setAttribute("class", "dynamicForm");
	var content;
	//console.log($xml);
	$($root.children()).each(function(i, node){
		//check value of each node, handle appropriately
		switch(node.nodeName){
			case "field":
				//get the form elements
				content = handleField(node);
				//append them
				form.append(content);
				break;
			case "fieldset":
				//get the form elements
				content = handleFieldset(node);
				//append them
				form.append(content);
				break;
			default:
				console.log("default");
		}
	})
	//console.log(form);
	//create the generate button for form
	var genButton = document.createElement("button");
	//make it do stuff
	genButton.innerHTML = "Generate Config";
	genButton.setAttribute("type", "button");
	genButton.setAttribute("class", "dynamicSubmit");
	genButton.setAttribute("onclick", "configGenerate()");
	//add to form
	form.append(genButton);
	
	//create wrapper for filename field
	var wrapper = document.createElement("div");
	//create label for filename field
	var fileNameLabel = document.createElement("label");
	fileNameLabel.innerHTML = "Filename to Save As"
	//create filename input field
	var fileNameField = document.createElement("input");
	fileNameField.setAttribute("type", "text");
	fileNameField.setAttribute("id", "fileNameField");
	fileNameField.setAttribute("name", "fileNameField");
	fileNameField.setAttribute("class", "dynamicFilename");
	//add label and input to wrapper
	wrapper.append(fileNameLabel);
	wrapper.append(fileNameField);
	
	document.getElementById("content").innerHTML = '';
	document.getElementById("content").append(wrapper);
	document.getElementById("content").append(form);
}

function handleField(node){
	//all nodes have a label
	//create node element
	var label = document.createElement('label');
	label.setAttribute("class", "dynamicLabel");
	//set the content of the label element to the text
		//within label xml
	label.innerHTML = $(node).find("label").text();
	//initialize the contentNode from xml
	var contentNode = $(node).children()[1];
	var content;
	var fieldFinal = document.createElement('div');
	//fields can be one of 4 types
	//call appropriate element creater based on type
	switch(contentNode.nodeName){
		case "inputText":
			//console.log(label+" is inputText");
			content = createTextField(contentNode);
			fieldFinal.setAttribute("class", "dynamicText");
			break;
		case "inputNumber":
			content = createNumberField(contentNode);
			fieldFinal.setAttribute("class", "dynamicNumber");
			break;
		case "select":
			content = createSelectField(contentNode);
			fieldFinal.setAttribute("class", "dynamicSelect");
			break;
		case "checkboxList":
			content = createCheckboxField(contentNode);
			fieldFinal.setAttribute("class", "dynamicCheckbox");
			label.setAttribute("class", "dynamicCheckboxLabel");
			break;
		default:
			console.log("Something went wrong.");
	}
	//append elements to the div
	fieldFinal.appendChild(label);
	fieldFinal.appendChild(content);
	//return the div
	//console.log(fieldFinal);
	return fieldFinal;
}

function handleFieldset(node){
	//console.log("fieldset");
	var content = document.createElement('fieldset');
	content.setAttribute("class", "dynamicFieldset");
	var legend = document.createElement('legend');
	//console.log($(node).attr('legend'));
	legend.innerHTML = $(node).attr('legend');
	content.append(legend);
	var tempContent;
	//legend.innerHTML = $(node).attr('legend');
	$(node).children().each(function(i, field){
		if(field.nodeName == 'fieldset'){
			tempContent = handleFieldset(field);
		} else {
			tempContent = handleField(field);
		}
		content.append(tempContent);
	})
	return content;
}
/*
function handleFieldsetRecursive(node){
	var content;
	if($(node).nodeName == "fieldset"){
		content = handleFieldset(node);
	} else {
		content = handleField(node);
	}
	return content;
}*/

function createTextField(node){
	var content = document.createElement('input');
	content.setAttribute("type", "text");
	content.setAttribute("id", $(node).attr('id'));
	content.setAttribute("name", $(node).attr('name'));
	return content;
}

function createNumberField(node){
	var content = document.createElement('input');
	content.setAttribute("type", "number");
	content.setAttribute("id", $(node).attr('id'));
	content.setAttribute("name", $(node).attr('name'));
	return content;
}

function createSelectField(node){
	var content = document.createElement('select');
	content.setAttribute("id", $(node).attr('id'));
	content.setAttribute("name", $(node).attr('name'));
	//for each child, create a select label option
	$(node).children().each(function(i, child){
		var option = document.createElement('option');
		option.setAttribute("value", $(child).attr('value'));
		option.innerHTML = $(child).text();
		content.append(option);
	})
	return content;
}

function createCheckboxField(node){
	//console.log(node.nodeName);
	var superContent = document.createElement('div');
	var content = document.createElement('ul');
	content.setAttribute("class", "dynamicCheckboxContent");
	//for each child, create a list item 
	//with a checkbox inside
	$(node).children().each(function(i, child){
		var wrapper = document.createElement('li');
		var checkboxItem = document.createElement('input');
		checkboxItem.setAttribute("type", "checkbox");
		checkboxItem.setAttribute("id", $(child).attr('id'));
		checkboxItem.setAttribute("name", $(node).attr('name'));
		checkboxItem.setAttribute("value", $(child).attr('value'));
		wrapper.append(checkboxItem);
		wrapper.append($(child).text());
		content.append(wrapper);
	})
	superContent.append(content);
	return superContent;
}

function configGenerate() {
	console.log("Something");
	//converts form fields into Java array
    var form = $("#content").children()[1];
	console.log("Help me please");
	console.log(form);
	var array = $(form).serializeArray();
    console.log(array);
    //perform string cleaning here
    var str = "";
    var fName = "";
	
	//iterate through each item in the array and append to str
    jQuery.each( array, function( i, field ) {
        //if name is same as last item, append to that line instead
        if(fName == field.name){
            str += field.value;
        } else {
            fName = field.name;
            str += "\n" + field.name + "=";
            str += field.value;
        }
    });
	
	//download the formatted string as a .cfg file
    download(document.getElementById("fileNameField").value + ".cfg", str);
}

function getFile() {
	$("#XMLFilename").click();
}


function updateText() {
	if(document.getElementById("XMLFilename").files[0]){
		var file = handleFileSelect('XMLFilename');
		var filename = file.name;
		console.log(filename);
		$("#placeText").html(filename);
	}
}

$("#XMLFilename").change(function() {
	if(document.getElementById("XMLFilename").files[0]){
		var file = handleFileSelect('XMLFilename');
		var filename = file.name;
		$("#placeText").html(filename);
	}
})

//download a string as contents of a file
function download(filename, text) {
    var element = document.createElement('a');
    element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(text));
    element.setAttribute('download', filename);

    element.style.display = 'none';
    document.body.appendChild(element);

    element.click();

    document.body.removeChild(element);
}