cwlVersion: v1.0
class: ExpressionTool
id: expression_preparerg
requirements:
  - class: InlineJavascriptRequirement
inputs:
  rg: {type: File, inputBinding: {loadContents: true}}
  sample: string
outputs:
  rg_str: string

expression:
  "${
      var arr = inputs.rg.contents.split('\\n')[0].split('\\t');
      for (var i=1; i<arr.length; i++){
        if (arr[i].startsWith('SM')){
          arr[i] = 'SM:' + inputs.sample;
        }
      }
      return {rg_str: arr.join('\\\\t')};
    }"
