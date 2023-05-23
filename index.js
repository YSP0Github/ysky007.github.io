let oSearchBar=document.querySelector('.searchBar');
let oIcon=document.querySelector('.icon');
let oClear=document.querySelector('.clear');
let oText=document.querySelector('input[type="text"]');

oIcon.addEventListener('click',() => {
    oSearchBar.classList.toggle('changeWidth')
});

oClear.addEventListener('click',() => {
    oText.value = '';
});