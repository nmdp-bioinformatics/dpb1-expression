import { Component, Input, OnInit } from '@angular/core';

@Component({
  selector: 'app-direction-match-status',
  templateUrl: './direction-match-status.component.html',
  styleUrls: ['./direction-match-status.component.scss']
})
export class DirectionMatchStatusComponent implements OnInit {
  @Input() data: string;

  constructor() { }

  ngOnInit() {
  }

}
