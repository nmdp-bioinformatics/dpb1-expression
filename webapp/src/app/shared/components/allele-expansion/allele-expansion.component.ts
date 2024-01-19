import { Component, Inject, OnInit } from '@angular/core';
import { MatDialogRef, MAT_DIALOG_DATA } from '@angular/material/dialog';
import { Allotype } from '@app/shared/models/allotype/allotype.model';

export interface DialogData {
  allotype : Allotype
}

@Component({
  selector: 'app-allele-expansion',
  templateUrl: './allele-expansion.component.html',
  styleUrls: ['./allele-expansion.component.scss']
})
export class AlleleExpansionComponent implements OnInit {

  constructor(
    public dialogRef: MatDialogRef<AlleleExpansionComponent>,
    @Inject(MAT_DIALOG_DATA) public data: DialogData,
  ) { }

  displayedColumns : string[] = ['typing', 'ciwd', 'exprLevel', 'experimental'];
  // 'exon3Motif',
  // 'utr3Motif', 

  ngOnInit() {
  }

}
